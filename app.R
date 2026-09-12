# =============================================================================
# GoldenEye v2.0 - Mussel Treatment-Tag Colonization Monitoring
# =============================================================================
# Created by L. McCormick (lmccormick@usbr.gov) with original concept and v1
# from Joe Hatfield (jhatfield@drreedinc.com)
# v2.0 developed 2026/09/11
#
# Single-file Shiny app organized into clearly separated module sections:
#   1. Pure calculation functions (calc_geometry, calc_coverage) - unit-testable
#   2. mod_metadata   - station/tag/treatment/material/depth/side inputs
#   3. mod_canvas     - JS canvas bridge (corners, ruler, cutout, paint)
#   4. mod_qa         - QA panel: invariant check + area deviation flag
#   5. mod_storage    - RDS-backed record read/write, CSV export
#   6. app.R (bottom) - ties modules together
#
# STORAGE: local RDS file (data/goldeneye_records.rds) + images saved to
# data/images/. This is intentionally simple for v1; swapping to a database
# later only requires rewriting mod_storage functions, since all other
# modules interact with storage only through save_record()/load_records().
#
# SCOPE NOTE: Blind mode is NOT implemented in this verison.
# The metadata module is structured so blind mode can be added later as a
# conditional branch without restructuring - see comment in mod_metadata_ui.
# =============================================================================

library(shiny)
library(shinyjs)
library(uuid)
library(jsonlite)
library(DT)
library(base64enc)

#' Normalizes a point list received from Shiny.setInputValue() JSON payload
#' into a guaranteed plain R list of list(x=numeric, y=numeric) elements.
#'
#' WHY THIS EXISTS: Shiny's default JSON deserialization (via jsonlite) applies
#' automatic simplification to arrays of objects - a JS array of {x,y} points
#' can arrive as a data.frame, a matrix, or a list depending on point count and
#' structure, rather than reliably as a list of named pairs. Every downstream
#' geometry function (polygon_area_px2, point_distance_px) expects a plain list
#' of list(x=,y=) elements. Without this normalization step, type mismatches
#' here were the root cause of the "is.character(txt) is not TRUE" errors -
#' malformed/unexpected shapes were propagating into validate()/need() calls
#' and other logic that assumed a specific list structure.
#'
#' @param raw the raw value from the payload (could be NULL, a data.frame,
#'   a list of lists, or a single list(x=,y=) if simplified to one point)
#' @return a plain list of list(x=,y=) elements, or an empty list if raw is
#'   NULL/empty/unparseable - NEVER errors, so callers can safely check length()
normalize_points <- function(raw) {
  if (is.null(raw)) return(list())
  
  # Case 1: already a plain list of list(x=,y=) elements (most common when
  # jsonlite does NOT simplify, e.g. with a mix of point counts across fields)
  if (is.list(raw) && !is.data.frame(raw)) {
    # Could be list(list(x=,y=), list(x=,y=), ...) OR a single list(x=,y=)
    if (!is.null(raw$x) && !is.null(raw$y) && is.numeric(raw$x) && length(raw$x) == 1) {
      # Single point object, e.g. ruler_point_1 - wrap it for consistency
      return(list(list(x = as.numeric(raw$x), y = as.numeric(raw$y))))
    }
    out <- lapply(raw, function(pt) {
      if (is.list(pt) && !is.null(pt$x) && !is.null(pt$y)) {
        list(x = as.numeric(pt$x), y = as.numeric(pt$y))
      } else {
        NULL
      }
    })
    out <- Filter(Negate(is.null), out)
    return(out)
  }
  
  # Case 2: jsonlite simplified the array into a data.frame with x, y columns
  if (is.data.frame(raw)) {
    if (!all(c("x", "y") %in% names(raw))) return(list())
    return(lapply(seq_len(nrow(raw)), function(i) {
      list(x = as.numeric(raw$x[i]), y = as.numeric(raw$y[i]))
    }))
  }
  
  # Case 3: jsonlite simplified into a matrix
  if (is.matrix(raw)) {
    if (!all(c("x", "y") %in% colnames(raw))) return(list())
    return(lapply(seq_len(nrow(raw)), function(i) {
      list(x = as.numeric(raw[i, "x"]), y = as.numeric(raw[i, "y"]))
    }))
  }
  
  # Unrecognized shape - fail safe rather than error
  list()
}

#' Normalizes a single point (e.g. ruler_point_1) into list(x=,y=) or NULL.
normalize_single_point <- function(raw) {
  if (is.null(raw)) return(NULL)
  pts <- normalize_points(raw)
  if (length(pts) == 0) return(NULL)
  pts[[1]]
}


# =============================================================================
# SECTION 1: PURE CALCULATION FUNCTIONS
# These have no Shiny dependencies and should be unit-tested independently
# (e.g., with testthat) exactly as specified in the handoff doc, section 9.
# =============================================================================

#' Shoelace formula for polygon area in pixel^2.
#' Works for any simple (non-self-intersecting) polygon: tag corners,
#' tab corners, or cutout polygon all use this same function.
#' @param points list of list(x=, y=) or data.frame with x,y columns
polygon_area_px2 <- function(points) {
  df <- as.data.frame(do.call(rbind, lapply(points, function(p) c(x = p$x, y = p$y))))
  n <- nrow(df)
  if (n < 3) return(0)
  x <- df$x
  y <- df$y
  x_next <- c(x[-1], x[1])
  y_next <- c(y[-1], y[1])
  area <- abs(sum(x * y_next - x_next * y)) / 2
  area
}

#' Euclidean distance between two points, in pixels.
point_distance_px <- function(p1, p2) {
  sqrt((p2$x - p1$x)^2 + (p2$y - p1$y)^2)
}

#' cm per rectified pixel, from ruler click points + known real-world distance.
#' @param ruler_point_1,ruler_point_2 list(x=, y=) in RECTIFIED coordinate space
#' @param ruler_known_distance_cm numeric, analyst-entered
cm_per_rectified_px <- function(ruler_point_1, ruler_point_2, ruler_known_distance_cm) {
  px_dist <- point_distance_px(ruler_point_1, ruler_point_2)
  if (px_dist <= 0) {
    stop("Ruler points are identical or invalid; cannot compute scale.")
  }
  ruler_known_distance_cm / px_dist
}

#' Compute tag area (cm^2), handling single-region (Aluminum) and
#' two-region composite (Plastic: body + tab), with optional cutout subtraction.
#'
#' @param body_corners list of 4 list(x=,y=) in rectified space
#' @param tab_corners list of 4 list(x=,y=) in rectified space, or NULL for Aluminum
#' @param cutout_polygon list of N list(x=,y=) in rectified space, or NULL if no cutout
#' @param cm_per_px numeric scale factor from cm_per_rectified_px()
#' @return list with body_area_cm2, tab_area_cm2, cutout_area_cm2, tag_area_cm2
compute_tag_area <- function(body_corners, tab_corners = NULL, cutout_polygon = NULL, cm_per_px) {
  body_area_cm2 <- polygon_area_px2(body_corners) * cm_per_px^2
  
  tab_area_cm2 <- if (!is.null(tab_corners) && length(tab_corners) >= 3) {
    polygon_area_px2(tab_corners) * cm_per_px^2
  } else {
    NA_real_
  }
  
  cutout_area_cm2 <- if (!is.null(cutout_polygon) && length(cutout_polygon) >= 3) {
    polygon_area_px2(cutout_polygon) * cm_per_px^2
  } else {
    0
  }
  
  tab_component <- if (is.na(tab_area_cm2)) 0 else tab_area_cm2
  tag_area_cm2 <- body_area_cm2 + tab_component - cutout_area_cm2
  
  list(
    body_area_cm2 = body_area_cm2,
    tab_area_cm2 = tab_area_cm2,
    cutout_area_cm2 = cutout_area_cm2,
    tag_area_cm2 = tag_area_cm2
  )
}

#' Compute coverage percentages from pixel-equivalent class totals.
#' Pixel-equivalents should already be SUMMED ACROSS REGIONS (body+tab) by
#' the caller before this function is invoked - this function does not know
#' about regions, only about total class areas, per spec section 4.1.1.
#'
#' @param mussel_px,bio_px,obscured_px pixel-equivalent totals (summed across regions)
#' @param total_tag_area_px total rectified tag area in px^2 (summed across regions,
#'   NOT net of cutout in pixel space - cutout is handled by excluding it from
#'   paintable area upstream, so obscured/mussel/bio pixel totals should already
#'   reflect a canvas that excludes the cutout region)
#' @return list of percentages + measurable_area_px + invariant check
compute_coverage_pct <- function(mussel_px, bio_px, obscured_px, total_tag_area_px,
                                 tolerance_pct = 0.05) {
  if (total_tag_area_px <= 0) stop("total_tag_area_px must be positive.")
  
  obscured_pct <- (obscured_px / total_tag_area_px) * 100
  measurable_area_px <- total_tag_area_px - obscured_px
  
  if (measurable_area_px <= 0) {
    warning("Measurable area is zero or negative - fully obscured tag?")
    mussel_pct <- NA_real_
    bio_pct <- NA_real_
    clean_pct <- NA_real_
  } else {
    mussel_pct <- (mussel_px / measurable_area_px) * 100
    bio_pct <- (bio_px / measurable_area_px) * 100
    unpainted_px <- measurable_area_px - mussel_px - bio_px
    clean_pct <- (unpainted_px / measurable_area_px) * 100
  }
  
  invariant_sum <- mussel_pct + bio_pct + clean_pct
  invariant_ok <- !is.na(invariant_sum) && abs(invariant_sum - 100) <= tolerance_pct
  
  list(
    obscured_pct = obscured_pct,
    mussel_pct = mussel_pct,
    biofouling_pct = bio_pct,
    clean_substrate_pct = clean_pct,
    measurable_area_px = measurable_area_px,
    invariant_sum = invariant_sum,
    invariant_ok = invariant_ok
  )
}

#' Estimated physical mussel area (cm^2).
compute_mussel_area_cm2 <- function(tag_area_cm2, obscured_pct, mussel_pct) {
  obscured_fraction <- obscured_pct / 100
  mussel_fraction <- mussel_pct / 100
  tag_area_cm2 * (1 - obscured_fraction) * mussel_fraction
}

#' QA flag: does computed tag_area_cm2 deviate more than threshold_pct from target?
#' Returns NA (not FALSE) if no target is configured, so callers can distinguish
#' "not flagged because within tolerance" from "not flagged because no target set".
check_area_qa_flag <- function(tag_area_cm2, target_cm2, threshold_pct = 0.05) {
  if (is.null(target_cm2) || is.na(target_cm2) || target_cm2 <= 0) {
    return(NA)
  }
  deviation <- abs(tag_area_cm2 - target_cm2) / target_cm2
  deviation > threshold_pct
}

#' Validate that all observations for a given tag_id share the same date_deployed.
#' @param records data.frame of existing records
#' @param tag_id character
#' @param new_date_deployed Date/character to check against existing records
#' @return list(ok = logical, existing_date = value or NA, message = character)
check_date_deployed_consistency <- function(records, tag_id, new_date_deployed) {
  existing <- records[records$tag_id == tag_id, ]
  if (nrow(existing) == 0) {
    return(list(ok = TRUE, existing_date = NA, message = "No prior records for this tag_id."))
  }
  existing_dates <- unique(as.character(existing$date_deployed))
  if (length(existing_dates) > 1) {
    return(list(
      ok = FALSE, existing_date = existing_dates,
      message = sprintf(
        "Existing records for tag_id '%s' already have INCONSISTENT date_deployed values: %s",
        tag_id, paste(existing_dates, collapse = ", ")
      )
    ))
  }
  if (as.character(new_date_deployed) != existing_dates[1]) {
    return(list(
      ok = FALSE, existing_date = existing_dates[1],
      message = sprintf(
        "date_deployed mismatch for tag_id '%s': existing records use %s, this entry uses %s.",
        tag_id, existing_dates[1], as.character(new_date_deployed)
      )
    ))
  }
  list(ok = TRUE, existing_date = existing_dates[1], message = "Consistent.")
}


# =============================================================================
# SECTION 2: mod_metadata - metadata entry UI/server
# =============================================================================

mod_metadata_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # NOTE ON BLIND MODE (deferred): when implemented, wrap the block below
    # in conditionalPanel(condition = "!input.blind_mode", ...) and add a
    # parallel minimal panel (filename + tag_id + side only) shown when
    # blind_mode is TRUE. Metadata resolution would then happen server-side
    # via a lookup table uploaded in a separate admin module, matched on
    # filename or tag_id, and joined onto the record after submission.
    textInput(ns("station_id"), "Station ID"),
    textInput(ns("tag_id"), "Tag ID (match photo notation)"),
    radioButtons(ns("side"), "Side", choices = c("A", "B"), inline = TRUE),
    radioButtons(ns("tag_material"), "Tag Material",
                 choices = c("Aluminum", "Plastic"), inline = TRUE),
    dateInput(ns("date_deployed"), "Date Deployed"),
    dateInput(ns("observation_date"), "Observation Date", value = Sys.Date()),
    timeInput <- textInput(ns("observation_time"), "Observation Time (HH:MM)",
                           value = format(Sys.time(), "%H:%M")),
    selectInput(ns("treatment"), "Treatment",
                choices = c("Control", "Barnacle Stop", "Nylon Stocking", "Other")),
    radioButtons(ns("reapply_schedule"), "Reapply Schedule",
                 choices = c("Monthly reapply", "No reapply (rinse only)")),
    numericInput(ns("depth"), "Depth (1 = shallowest, 12 = deepest)",
                 value = 1, min = 1, max = 12, step = 1),
    textInput(ns("analyst"), "Analyst Name")
  )
}

mod_metadata_server <- function(id, existing_records) {
  moduleServer(id, function(input, output, session) {
    
    # Live consistency check against existing records as tag_id/date change.
    date_check <- reactive({
      req(input$tag_id, input$date_deployed)
      check_date_deployed_consistency(existing_records(), input$tag_id, input$date_deployed)
    })
    
    output$date_consistency_msg <- renderText({
      dc <- date_check()
      if (!dc$ok) paste("WARNING:", dc$message) else ""
    })
    
    # Bundles all metadata into one list for downstream use.
    metadata <- reactive({
      req(input$station_id, input$tag_id, input$date_deployed)
      list(
        station_id = input$station_id,
        tag_id = input$tag_id,
        side = input$side,
        tag_material = input$tag_material,
        date_deployed = as.character(input$date_deployed),
        observation_datetime = paste(input$observation_date, input$observation_time),
        treatment = input$treatment,
        reapply_schedule = input$reapply_schedule,
        depth = input$depth,
        analyst = input$analyst
      )
    })
    
    list(
      metadata = metadata,
      date_check = date_check
    )
  })
}


# =============================================================================
# SECTION 3: mod_canvas - JS canvas bridge
# =============================================================================
# This is the most substantial JS/R integration point. The JS side handles
# all interactive drawing (corners, ruler, cutout, paint/erase) and sends a
# structured payload via Shiny.setInputValue(). This implementation uses a
# hidden offscreen mask canvas with flat, fully opaque per-class colors and
# getImageData() to count actual unique painted pixels (native Canvas2D API,
# no external libraries) - this avoids double-counting overlapping strokes,
# which a geometric "sum of brush circles" approach cannot do correctly.
# =============================================================================

mod_canvas_ui <- function(id) {
  ns <- NS(id)
  tagList(
    useShinyjs(),
    fileInput(ns("image_upload"), "Upload Tag Photo", accept = c("image/jpeg", "image/png")),
    
    tags$div(
      style = "margin-bottom: 8px;",
      tags$strong("Step: "),
      tags$span(id = ns("current_step_label"), "Upload an image to begin")
    ),
    
    tags$div(
      style = "display:flex; gap:8px; flex-wrap: wrap; margin-bottom:8px;",
      actionButton(ns("mode_body"), "1. Select Body Corners"),
      actionButton(ns("mode_tab"), "2. Select Tab Corners (Plastic only)"),
      actionButton(ns("mode_ruler"), "3. Mark Ruler"),
      actionButton(ns("mode_cutout"), "4. Trace Cutout (optional)"),
      actionButton(ns("finish_cutout"), "Finish Cutout Shape"),
      actionButton(ns("mode_paint_mussel"), "5. Paint: Mussel",
                   style = "background-color: rgb(230,159,0); color: white; border: none;"),
      actionButton(ns("mode_paint_bio"), "Paint: Biofouling",
                   style = "background-color: rgb(0,114,178); color: white; border: none;"),
      actionButton(ns("mode_paint_obscured"), "Paint: Obscured",
                   style = "background-color: rgb(204,121,167); color: white; border: none;"),
      actionButton(ns("mode_erase"), "Erase")
    ),
    
    conditionalPanel(
      condition = "input.mode_ruler > 0",
      ns = ns,
      tags$div(
        style = "display:flex; align-items:center; gap:8px;",
        numericInput(ns("ruler_known_distance_value"), "Ruler distance represents:",
                     value = 10, min = 0.01, step = 0.1),
        radioButtons(ns("ruler_unit"), NULL, choices = c("cm", "in"),
                     selected = "cm", inline = TRUE)
      ),
      helpText("Enter the value in whichever unit is convenient - it is converted to ",
               "centimeters automatically before any calculation. All stored and displayed ",
               "areas remain in metric (cm\u00b2).")
    ),
    
    tags$div(
      style = "display:flex; align-items:center; gap:12px; margin-bottom:8px;",
      tags$label("Mask opacity:", style = "margin-bottom:0;"),
      tags$input(id = ns("mask_opacity"), type = "range", min = "0", max = "100",
                 value = "50", style = "width:200px;"),
      tags$span(id = ns("mask_opacity_label"), "50%"),
      tags$button(id = ns("toggle_hide_mask"), type = "button",
                  class = "btn btn-outline-secondary",
                  "Hold to Hide Mask (check your work)")
    ),
    
    tags$div(
      tags$canvas(id = ns("main_canvas"), width = 1200, height = 800,
                  style = "border:1px solid #999; max-width:100%; cursor: crosshair;"),
      # Hidden mask canvas - never displayed, used only for pixel-accurate class counting.
      tags$canvas(id = ns("mask_canvas"), width = 1200, height = 800,
                  style = "display:none;")
    ),
    
    actionButton(ns("reset_canvas"), "Reset / Start Over", class = "btn-danger"),
    actionButton(ns("send_to_r"), "Calculate Coverage", class = "btn-primary"),
    
    tags$script(HTML(mod_canvas_js(ns("main_canvas"), ns("mask_canvas"), ns(""))))
  )
}

#' Generates the canvas JS as a plain string. Uses simple token replacement
#' (gsub with fixed=TRUE) rather than sprintf(), because the JS body contains
#' literal '%' characters (e.g. opacity label text, CSS-adjacent values) that
#' sprintf() would otherwise misinterpret as format specifiers and error on.
#' @param canvas_id full namespaced DOM id of the visible canvas
#' @param mask_canvas_id full namespaced DOM id of the hidden mask canvas
#' @param ns_prefix namespace prefix used for all other input/button ids,
#'   e.g. "canvas-" (module id "canvas" produces ids like "canvas-mode_body")
mod_canvas_js <- function(canvas_id, mask_canvas_id, ns_prefix) {
  js <- "
(function() {
  const canvasId = '__CANVAS_ID__';
  const maskCanvasId = '__MASK_CANVAS_ID__';
  const nsPrefix = '__NS_PREFIX__'; // e.g. 'canvas-'

  let canvas, ctx, maskCanvas, maskCtx;
  let img = null;
  let mode = null; // 'body' | 'tab' | 'ruler' | 'cutout' | 'paint_mussel' | 'paint_bio' | 'paint_obscured' | 'erase'
  let bodyCorners = [];
  let tabCorners = [];
  let rulerPoints = [];
  let cutoutPoints = [];
  let cutoutFinished = false;
  let maskOpacity = 0.5; // live, analyst-adjustable; VISUAL ONLY - never affects counting
  let maskHidden = false; // temporary full-hide toggle; VISUAL ONLY - never affects counting
  const BRUSH_RADIUS = 12;

  // Flat, fully opaque colors used ONLY on the hidden mask canvas for counting.
  // Based on the Okabe-Ito colorblind-safe palette: orange and blue sit far
  // apart on every major axis of color vision deficiency (unlike orange/green,
  // which compress together under red-green color blindness, the most common
  // form). Reddish-purple for obscured is likewise off the red-green axis.
  const MASK_COLORS = {
    mussel:   'rgba(230,159,0,255)',   // Okabe-Ito orange
    bio:      'rgba(0,114,178,255)',   // Okabe-Ito blue
    obscured: 'rgba(204,121,167,255)'  // Okabe-Ito reddish purple
  };
  const MASK_RGBA = {
    mussel:   [230,159,0,255],
    bio:      [0,114,178,255],
    obscured: [204,121,167,255]
  };

  function id(suffix) { return nsPrefix + suffix; }

  function initCanvas() {
    canvas = document.getElementById(canvasId);
    ctx = canvas.getContext('2d');
    maskCanvas = document.getElementById(maskCanvasId);
    maskCtx = maskCanvas.getContext('2d', { willReadFrequently: true });
    maskCtx.clearRect(0, 0, maskCanvas.width, maskCanvas.height);
  }

  function redraw() {
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    if (img) ctx.drawImage(img, 0, 0, canvas.width, canvas.height);

    // Composite the mask on top for VISUAL feedback only, at analyst-controlled
    // opacity (or fully hidden via the hide toggle). Counting always reads from
    // maskCanvas directly via getImageData on its own opaque pixel data, so
    // nothing about this display step can ever bias the actual measurements -
    // this is purely so the analyst can see the underlying photo while painting.
    if (!maskHidden) {
      ctx.save();
      ctx.globalAlpha = maskOpacity;
      ctx.drawImage(maskCanvas, 0, 0);
      ctx.restore();
    }

    drawPolygon(bodyCorners, 'lime', true);
    drawPolygon(tabCorners, 'cyan', true);
    drawPolygon(cutoutPoints, 'red', cutoutFinished);

    if (rulerPoints.length > 0) {
      ctx.fillStyle = 'yellow';
      rulerPoints.forEach(p => {
        ctx.beginPath();
        ctx.arc(p.x, p.y, 5, 0, 2*Math.PI);
        ctx.fill();
      });
      if (rulerPoints.length === 2) {
        ctx.strokeStyle = 'yellow';
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.moveTo(rulerPoints[0].x, rulerPoints[0].y);
        ctx.lineTo(rulerPoints[1].x, rulerPoints[1].y);
        ctx.stroke();
      }
    }
  }

  function drawPolygon(points, color, closeIt) {
    if (points.length === 0) return;
    ctx.strokeStyle = color;
    ctx.fillStyle = color;
    ctx.lineWidth = 2;
    points.forEach(p => {
      ctx.beginPath();
      ctx.arc(p.x, p.y, 4, 0, 2*Math.PI);
      ctx.fill();
    });
    if (points.length > 1) {
      ctx.beginPath();
      ctx.moveTo(points[0].x, points[0].y);
      for (let i = 1; i < points.length; i++) ctx.lineTo(points[i].x, points[i].y);
      if (closeIt && points.length >= 3) ctx.closePath();
      ctx.stroke();
    }
  }

  function pointInPolygon(pt, poly) {
    if (poly.length < 3) return false;
    let inside = false;
    for (let i = 0, j = poly.length - 1; i < poly.length; j = i++) {
      const xi = poly[i].x, yi = poly[i].y;
      const xj = poly[j].x, yj = poly[j].y;
      const intersect = ((yi > pt.y) !== (yj > pt.y)) &&
        (pt.x < (xj - xi) * (pt.y - yi) / (yj - yi) + xi);
      if (intersect) inside = !inside;
    }
    return inside;
  }

  function getCanvasCoords(evt) {
    const rect = canvas.getBoundingClientRect();
    const scaleX = canvas.width / rect.width;
    const scaleY = canvas.height / rect.height;
    return {
      x: (evt.clientX - rect.left) * scaleX,
      y: (evt.clientY - rect.top) * scaleY
    };
  }

  function handleClick(evt) {
    const pt = getCanvasCoords(evt);
    if (mode === 'body' && bodyCorners.length < 4) {
      bodyCorners.push(pt);
    } else if (mode === 'tab' && tabCorners.length < 4) {
      tabCorners.push(pt);
    } else if (mode === 'ruler' && rulerPoints.length < 2) {
      rulerPoints.push(pt);
    } else if (mode === 'cutout' && !cutoutFinished) {
      cutoutPoints.push(pt);
    }
    redraw();
  }

  let isPainting = false;
  function handlePaintStart(evt) {
    if (!['paint_mussel','paint_bio','paint_obscured','erase'].includes(mode)) return;
    isPainting = true;
    handlePaintMove(evt);
  }

  function handlePaintMove(evt) {
    if (!isPainting) return;
    const pt = getCanvasCoords(evt);

    // Never allow painting inside a finished cutout region - it is not tag surface.
    if (cutoutFinished && pointInPolygon(pt, cutoutPoints)) return;

    if (mode === 'erase') {
      maskCtx.save();
      maskCtx.globalCompositeOperation = 'destination-out';
      maskCtx.beginPath();
      maskCtx.arc(pt.x, pt.y, BRUSH_RADIUS, 0, 2*Math.PI);
      maskCtx.fill();
      maskCtx.restore();
      redraw();
      return;
    }

    const cls = mode.replace('paint_', '');
    maskCtx.fillStyle = MASK_COLORS[cls];
    maskCtx.beginPath();
    maskCtx.arc(pt.x, pt.y, BRUSH_RADIUS, 0, 2*Math.PI);
    maskCtx.fill();
    redraw();
  }

  function handlePaintEnd() { isPainting = false; }

  function setMode(m) { mode = m; }

  function loadImageFile(file) {
    const reader = new FileReader();
    reader.onload = function(e) {
      img = new Image();
      img.onload = function() { redraw(); };
      img.src = e.target.result;
    };
    reader.readAsDataURL(file);
  }

  // Pixel-accurate class counting via getImageData - counts ACTUAL unique
  // painted pixels, so overlapping strokes of the same or different classes
  // are never double-counted. Optionally restricted to a region polygon
  // (Body vs Tab) by checking pointInPolygon per pixel.
  function countClassPixels(cls, regionPolygon) {
    const w = maskCanvas.width, h = maskCanvas.height;
    const imgData = maskCtx.getImageData(0, 0, w, h).data;
    const target = MASK_RGBA[cls];
    let count = 0;
    for (let y = 0; y < h; y++) {
      for (let x = 0; x < w; x++) {
        const idx = (y * w + x) * 4;
        if (imgData[idx] === target[0] && imgData[idx+1] === target[1] &&
            imgData[idx+2] === target[2] && imgData[idx+3] === target[3]) {
          if (!regionPolygon || pointInPolygon({x:x,y:y}, regionPolygon)) {
            count++;
          }
        }
      }
    }
    return count;
  }

  function assemblePayload() {
    const tabActive = tabCorners.length === 4;

    // Read the ruler distance value + unit and convert to cm BEFORE sending
    // to R, so every downstream calculation only ever deals with centimeters
    // regardless of what unit the analyst found convenient to type in.
    const rulerValueInput = document.getElementById(id('ruler_known_distance_value'));
    const rulerUnitRadios = document.querySelectorAll('input[name=\"' + id('ruler_unit') + '\"]');
    let rulerUnit = 'cm';
    rulerUnitRadios.forEach(r => { if (r.checked) rulerUnit = r.value; });
    const rulerRawValue = parseFloat(rulerValueInput?.value || '10');
    const INCHES_TO_CM = 2.54;
    const rulerDistanceCm = rulerUnit === 'in' ? rulerRawValue * INCHES_TO_CM : rulerRawValue;

    return {
      body_corners: bodyCorners,
      tab_corners: tabActive ? tabCorners : null,
      ruler_point_1: rulerPoints[0] || null,
      ruler_point_2: rulerPoints[1] || null,
      ruler_known_distance_cm: rulerDistanceCm,
      ruler_input_value_raw: rulerRawValue,
      ruler_input_unit_raw: rulerUnit,
      cutout_polygon: (cutoutFinished && cutoutPoints.length >= 3) ? cutoutPoints : null,
      // Pixel-accurate counts, split by region if a tab exists, else all counted as Body.
      mussel_px_body: countClassPixels('mussel', tabActive ? bodyCorners : null),
      mussel_px_tab: tabActive ? countClassPixels('mussel', tabCorners) : 0,
      bio_px_body: countClassPixels('bio', tabActive ? bodyCorners : null),
      bio_px_tab: tabActive ? countClassPixels('bio', tabCorners) : 0,
      obscured_px_body: countClassPixels('obscured', tabActive ? bodyCorners : null),
      obscured_px_tab: tabActive ? countClassPixels('obscured', tabCorners) : 0,
      source_width: img ? img.naturalWidth : null,
      source_height: img ? img.naturalHeight : null,
      canvas_width: canvas.width,
      canvas_height: canvas.height,
      marked_image_data_url: canvas.toDataURL('image/png')
    };
  }

  function resetAll() {
    bodyCorners = []; tabCorners = []; rulerPoints = []; cutoutPoints = [];
    cutoutFinished = false;
    maskCtx.clearRect(0, 0, maskCanvas.width, maskCanvas.height);
    redraw();
  }

  document.addEventListener('DOMContentLoaded', function() {
    initCanvas();
    canvas.addEventListener('click', handleClick);
    canvas.addEventListener('mousedown', handlePaintStart);
    canvas.addEventListener('mousemove', handlePaintMove);
    window.addEventListener('mouseup', handlePaintEnd);

    document.getElementById(id('image_upload')).addEventListener('change', function(e) {
      if (e.target.files[0]) loadImageFile(e.target.files[0]);
    });
    document.getElementById(id('mode_body')).addEventListener('click', () => setMode('body'));
    document.getElementById(id('mode_tab')).addEventListener('click', () => setMode('tab'));
    document.getElementById(id('mode_ruler')).addEventListener('click', () => setMode('ruler'));
    document.getElementById(id('mode_cutout')).addEventListener('click', () => setMode('cutout'));
    document.getElementById(id('finish_cutout')).addEventListener('click', function() {
      if (cutoutPoints.length >= 3) { cutoutFinished = true; redraw(); }
    });
    document.getElementById(id('mode_paint_mussel')).addEventListener('click', () => setMode('paint_mussel'));
    document.getElementById(id('mode_paint_bio')).addEventListener('click', () => setMode('paint_bio'));
    document.getElementById(id('mode_paint_obscured')).addEventListener('click', () => setMode('paint_obscured'));
    document.getElementById(id('mode_erase')).addEventListener('click', () => setMode('erase'));

    document.getElementById(id('reset_canvas')).addEventListener('click', resetAll);

    // Live opacity slider - purely visual, updates on every drag frame.
    const opacitySlider = document.getElementById(id('mask_opacity'));
    const opacityLabel = document.getElementById(id('mask_opacity_label'));
    if (opacitySlider) {
      opacitySlider.addEventListener('input', function(e) {
        maskOpacity = parseInt(e.target.value, 10) / 100;
        if (opacityLabel) opacityLabel.textContent = e.target.value + '%';
        redraw();
      });
    }

    // Hide-mask toggle: press-and-hold to peek at the raw photo, release to
    // resume painting view. Works with both mouse (mousedown/mouseup) and
    // touch, and also supports a plain click as a fallback toggle for
    // accessibility (some analysts may find press-and-hold awkward).
    const hideBtn = document.getElementById(id('toggle_hide_mask'));
    if (hideBtn) {
      let heldViaMouse = false;
      hideBtn.addEventListener('mousedown', function() {
        heldViaMouse = true;
        maskHidden = true;
        redraw();
      });
      hideBtn.addEventListener('mouseup', function() {
        if (heldViaMouse) { maskHidden = false; heldViaMouse = false; redraw(); }
      });
      hideBtn.addEventListener('mouseleave', function() {
        if (heldViaMouse) { maskHidden = false; heldViaMouse = false; redraw(); }
      });
    }

    document.getElementById(id('send_to_r')).addEventListener('click', function() {
      const payload = assemblePayload();
      Shiny.setInputValue(nsPrefix + 'canvas_payload', payload, {priority: 'event'});
    });
  });
})();
"
  js <- gsub("__CANVAS_ID__", canvas_id, js, fixed = TRUE)
  js <- gsub("__MASK_CANVAS_ID__", mask_canvas_id, js, fixed = TRUE)
  js <- gsub("__NS_PREFIX__", ns_prefix, js, fixed = TRUE)
  js
}

mod_canvas_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    payload <- reactive({
      req(input$canvas_payload)
      input$canvas_payload
    })
    list(payload = payload)
  })
}


# =============================================================================
# SECTION 4: mod_qa - QA panel
# =============================================================================

mod_qa_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$h4("QA Status"),
    tableOutput(ns("qa_table")),
    numericInput(ns("target_aluminum_cm2"), "Target area - Aluminum (cm²)", value = NA, min = 0),
    numericInput(ns("target_plastic_cm2"), "Target area - Plastic (cm²)", value = NA, min = 0),
    helpText("Targets are manually configured, not computed from data (per spec section 7). ",
             "Leave blank to disable the QA flag until a target is set.")
  )
}

mod_qa_server <- function(id, tag_material, tag_area_cm2, invariant_ok, invariant_sum) {
  moduleServer(id, function(input, output, session) {
    
    area_flag <- reactive({
      mat <- tag_material()
      area <- tag_area_cm2()
      # Guard explicitly rather than relying on req() alone: if either input
      # isn't a valid length-1 value yet, return NA (not NULL) so downstream
      # code always has a consistent, checkable type to work with.
      if (is.null(mat) || length(mat) != 1 || is.na(mat) ||
          is.null(area) || length(area) != 1 || is.na(area)) {
        return(NA)
      }
      target <- if (mat == "Aluminum") input$target_aluminum_cm2 else input$target_plastic_cm2
      flag <- check_area_qa_flag(area, target, threshold_pct = 0.05)
      if (is.null(flag) || length(flag) != 1) NA else flag
    })
    
    output$qa_table <- renderTable({
      area <- tag_area_cm2()
      inv_sum <- invariant_sum()
      inv_ok <- invariant_ok()
      flag <- area_flag()
      
      # Build every Value entry as an explicit, guaranteed-length-1 character
      # string. This is defensive on purpose: renderTable's HTML generation
      # requires a character vector, and any NULL/NA/logical(0) slipping
      # through here (e.g. during a reactive's transient invalid state right
      # after a button click) previously caused an opaque
      # "is.character(txt) is not TRUE" error with no useful traceback.
      area_str <- if (is.null(area) || length(area) != 1 || is.na(area)) {
        "Pending calculation"
      } else {
        sprintf("%.2f", area)
      }
      
      flag_str <- if (is.null(flag) || length(flag) != 1 || is.na(flag)) {
        "No target set"
      } else if (isTRUE(flag)) {
        "FLAGGED (>5% deviation)"
      } else {
        "OK"
      }
      
      inv_sum_str <- if (is.null(inv_sum) || length(inv_sum) != 1 || is.na(inv_sum)) {
        "Pending calculation"
      } else {
        sprintf("%.3f", inv_sum)
      }
      
      inv_ok_str <- if (is.null(inv_ok) || length(inv_ok) != 1 || is.na(inv_ok)) {
        "Pending calculation"
      } else if (isTRUE(inv_ok)) {
        "OK"
      } else {
        "FAILED - investigate mask data"
      }
      
      data.frame(
        Check = c("Tag area (cm\u00b2)", "Area QA flag", "Coverage invariant sum, pct", "Invariant check"),
        Value = c(area_str, flag_str, inv_sum_str, inv_ok_str),
        stringsAsFactors = FALSE
      )
    }, sanitize.text.function = function(x) x)
    
    list(area_flag = area_flag)
  })
}


# =============================================================================
# SECTION 5: mod_storage - local RDS-backed storage
# =============================================================================

STORAGE_DIR <- "data"
RECORDS_PATH <- file.path(STORAGE_DIR, "goldeneye_records.rds")
IMAGES_DIR <- file.path(STORAGE_DIR, "images")

ensure_storage_dirs <- function() {
  if (!dir.exists(STORAGE_DIR)) dir.create(STORAGE_DIR, recursive = TRUE)
  if (!dir.exists(IMAGES_DIR)) dir.create(IMAGES_DIR, recursive = TRUE)
}

load_records <- function() {
  ensure_storage_dirs()
  if (!file.exists(RECORDS_PATH)) {
    return(data.frame(
      analysis_id = character(), tag_id = character(), side = character(),
      station_id = character(), tag_material = character(),
      treatment = character(), reapply_schedule = character(),
      depth = numeric(), date_deployed = character(), observation_datetime = character(),
      body_area_cm2 = numeric(), tab_area_cm2 = numeric(), cutout_area_cm2 = numeric(),
      tag_area_cm2 = numeric(), tag_area_qa_flag = logical(),
      mussel_cover_pct = numeric(), biofouling_cover_pct = numeric(),
      clean_substrate_pct = numeric(), obscured_pct = numeric(),
      mussel_area_cm2 = numeric(), blind_mode_used = logical(),
      analyst = character(), review_status = character(),
      created_at = character(), software_version = character(),
      ruler_input_value_raw = numeric(), ruler_input_unit_raw = character(),
      marked_image_path = character(),
      stringsAsFactors = FALSE
    ))
  }
  readRDS(RECORDS_PATH)
}

save_record <- function(record_row) {
  ensure_storage_dirs()
  existing <- load_records()
  updated <- rbind(existing, as.data.frame(record_row, stringsAsFactors = FALSE))
  saveRDS(updated, RECORDS_PATH)
  invisible(updated)
}

save_marked_image <- function(data_url, analysis_id) {
  ensure_storage_dirs()
  # data_url format: "data:image/png;base64,...."
  # Using base64enc::base64decode explicitly (not jsonlite) since that is the
  # function verified to exist with this exact name/signature - avoids relying
  # on an unverified cross-package function name for something this important.
  b64 <- sub("^data:image/png;base64,", "", data_url)
  raw <- base64enc::base64decode(b64)
  path <- file.path(IMAGES_DIR, paste0(analysis_id, ".png"))
  writeBin(as.raw(raw), path)
  path
}

export_csv <- function(path = "goldeneye_export.csv") {
  records <- load_records()
  write.csv(records, path, row.names = FALSE)
  path
}


# =============================================================================
# SECTION 6: app.R - main app tying modules together
# =============================================================================

ui <- fluidPage(
  useShinyjs(),
  titlePanel("GoldenEye v2.0 - Mussel Colonization Analysis"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      tags$h4("Metadata"),
      mod_metadata_ui("meta"),
      textOutput("date_consistency_warning"),
      tags$hr(),
      mod_qa_ui("qa")
    ),
    mainPanel(
      width = 9,
      mod_canvas_ui("canvas"),
      tags$hr(),
      tags$h4("Results Dashboard"),
      tableOutput("results_table"),
      actionButton("approve_record", "Approve & Save Record", class = "btn-success"),
      tags$hr(),
      tags$h4("Saved Records (this session)"),
      DTOutput("records_table"),
      downloadButton("download_csv", "Export CSV")
    )
  )
)

server <- function(input, output, session) {
  
  # Defined early, before any reactive that uses it - R does not hoist
  # function definitions the way some languages do, so this must come
  # before calc_results() below, which relies on it.
  `%||%` <- function(a, b) if (is.null(a)) b else a
  
  records_reactive <- reactiveVal(load_records())
  
  meta <- mod_metadata_server("meta", existing_records = records_reactive)
  canvas <- mod_canvas_server("canvas")
  
  output$date_consistency_warning <- renderText({
    dc <- meta$date_check()
    if (!is.null(dc) && !dc$ok) paste("\u26A0", dc$message) else ""
  })
  
  # --- Core calculation pipeline: canvas payload -> geometry -> coverage ---
  calc_results <- reactive({
    p <- canvas$payload()
    md <- meta$metadata()
    req(p, md)
    
    # Normalize every point field from the raw JSON payload BEFORE any
    # length()/validate() checks run. This is the fix for the
    # "is.character(txt) is not TRUE" error: raw JSON arrays of {x,y} points
    # can be deserialized by jsonlite as data.frames, matrices, or lists
    # depending on shape, and validate()/need() calls downstream were
    # receiving inconsistent types as a result. normalize_points() always
    # returns a plain list, so length() and every geometry function below
    # gets a guaranteed, consistent shape to work with.
    body_corners <- normalize_points(p$body_corners)
    tab_corners_raw <- normalize_points(p$tab_corners)
    tab_corners <- if (length(tab_corners_raw) == 4) tab_corners_raw else NULL
    cutout_polygon_raw <- normalize_points(p$cutout_polygon)
    cutout_polygon <- if (length(cutout_polygon_raw) >= 3) cutout_polygon_raw else NULL
    ruler_point_1 <- normalize_single_point(p$ruler_point_1)
    ruler_point_2 <- normalize_single_point(p$ruler_point_2)
    
    # Explicit, single validation block using plain if/stop with guaranteed
    # character messages, rather than chained validate()/need() calls - this
    # is more verbose but each message is constructed directly as a literal
    # string, so there is no path by which a non-character value could reach
    # Shiny's internal rendering and trigger the is.character(txt) assertion.
    validation_errors <- character(0)
    if (length(body_corners) != 4) {
      validation_errors <- c(validation_errors, "Please select all 4 body corners before calculating.")
    }
    if (is.null(ruler_point_1) || is.null(ruler_point_2)) {
      validation_errors <- c(validation_errors, "Please mark both ruler points before calculating.")
    }
    if (identical(md$tag_material, "Plastic") && is.null(tab_corners)) {
      validation_errors <- c(validation_errors, "Plastic tags require 4 tab corners in addition to body corners.")
    }
    if (length(validation_errors) > 0) {
      validate(need(FALSE, paste(validation_errors, collapse = " ")))
    }
    
    cm_per_px <- cm_per_rectified_px(ruler_point_1, ruler_point_2, p$ruler_known_distance_cm)
    
    area_result <- compute_tag_area(
      body_corners = body_corners,
      tab_corners = tab_corners,
      cutout_polygon = cutout_polygon,
      cm_per_px = cm_per_px
    )
    
    # Sum pixel-equivalents across regions BEFORE computing percentages,
    # per spec section 4.1.1 - single combined invariant check, not two.
    mussel_px <- (p$mussel_px_body %||% 0) + (p$mussel_px_tab %||% 0)
    bio_px <- (p$bio_px_body %||% 0) + (p$bio_px_tab %||% 0)
    obscured_px <- (p$obscured_px_body %||% 0) + (p$obscured_px_tab %||% 0)
    
    total_tag_area_px <- polygon_area_px2(body_corners) +
      (if (!is.null(tab_corners)) polygon_area_px2(tab_corners) else 0)
    
    coverage <- compute_coverage_pct(mussel_px, bio_px, obscured_px, total_tag_area_px)
    
    mussel_area_cm2 <- compute_mussel_area_cm2(
      area_result$tag_area_cm2, coverage$obscured_pct, coverage$mussel_pct
    )
    
    c(area_result, coverage, list(mussel_area_cm2 = mussel_area_cm2))
  })
  
  qa <- mod_qa_server(
    "qa",
    tag_material = reactive(meta$metadata()$tag_material),
    tag_area_cm2 = reactive(calc_results()$tag_area_cm2),
    invariant_ok = reactive(calc_results()$invariant_ok),
    invariant_sum = reactive(calc_results()$invariant_sum)
  )
  
  output$results_table <- renderTable({
    r <- calc_results()
    data.frame(
      Metric = c("Tag Area (cm²)", "Body Area (cm²)", "Tab Area (cm²)", "Cutout Area (cm²)",
                 "Mussel Cover (%)", "Biofouling Cover (%)", "Clean Substrate (%)",
                 "Obscured (%)", "Est. Mussel Area (cm²)"),
      Value = c(
        sprintf("%.2f", r$tag_area_cm2),
        sprintf("%.2f", r$body_area_cm2),
        ifelse(is.na(r$tab_area_cm2), "N/A", sprintf("%.2f", r$tab_area_cm2)),
        sprintf("%.2f", r$cutout_area_cm2),
        sprintf("%.2f", r$mussel_pct),
        sprintf("%.2f", r$biofouling_pct),
        sprintf("%.2f", r$clean_substrate_pct),
        sprintf("%.2f", r$obscured_pct),
        sprintf("%.2f", r$mussel_area_cm2)
      )
    )
  })
  
  observeEvent(input$approve_record, {
    # Entire save operation is wrapped in tryCatch so that ANY failure - a
    # missing package function, a bad path, malformed payload data, etc. -
    # surfaces as an on-screen modal instead of ever being able to propagate
    # up as an uncaught error. An uncaught error inside an observeEvent can,
    # in some Shiny configurations, terminate the session/websocket rather
    # than displaying Shiny's usual red error box, which is almost certainly
    # what caused the app window to close with no visible error previously.
    tryCatch({
      r <- calc_results()
      md <- meta$metadata()
      dc <- meta$date_check()
      
      if (!dc$ok) {
        showModal(modalDialog(
          title = "Date Deployed Inconsistency",
          dc$message,
          "Please correct date_deployed before saving, or confirm this is intentional.",
          footer = modalButton("Close")
        ))
        return(invisible(NULL))
      }
      
      analysis_id <- uuid::UUIDgenerate()
      p <- canvas$payload()
      
      marked_image_path <- tryCatch(
        save_marked_image(p$marked_image_data_url, analysis_id),
        error = function(e) {
          showNotification(
            paste("Warning: could not save marked image -", conditionMessage(e)),
            type = "warning", duration = NULL
          )
          NA_character_
        }
      )
      
      new_record <- list(
        analysis_id = analysis_id,
        tag_id = md$tag_id,
        side = md$side,
        station_id = md$station_id,
        tag_material = md$tag_material,
        treatment = md$treatment,
        reapply_schedule = md$reapply_schedule,
        depth = md$depth,
        date_deployed = md$date_deployed,
        observation_datetime = md$observation_datetime,
        body_area_cm2 = r$body_area_cm2,
        tab_area_cm2 = ifelse(is.na(r$tab_area_cm2), NA_real_, r$tab_area_cm2),
        cutout_area_cm2 = r$cutout_area_cm2,
        tag_area_cm2 = r$tag_area_cm2,
        tag_area_qa_flag = isTRUE(qa$area_flag()),
        mussel_cover_pct = r$mussel_pct,
        biofouling_cover_pct = r$biofouling_pct,
        clean_substrate_pct = r$clean_substrate_pct,
        obscured_pct = r$obscured_pct,
        mussel_area_cm2 = r$mussel_area_cm2,
        blind_mode_used = FALSE,
        analyst = md$analyst,
        review_status = "Draft",
        created_at = as.character(Sys.time()),
        software_version = "GoldenEye v2.0",
        marked_image_path = if (is.na(marked_image_path)) NA_character_ else marked_image_path
      )
      
      save_record(new_record)
      records_reactive(load_records())
      
      showNotification("Record saved.", type = "message")
      
    }, error = function(e) {
      # traceback() here writes the full call stack to the R console (not the
      # browser) so we can see exactly which function call produced the error,
      # rather than only seeing the final, unhelpful "is.character(txt)" message.
      message("=== approve_record error ===")
      message(conditionMessage(e))
      print(sys.calls())
      
      showModal(modalDialog(
        title = "Error Saving Record",
        tags$p("The record could not be saved. Your data has NOT been lost from this session -"),
        tags$p("the form and canvas still hold your work. Full details were also printed to the R console."),
        tags$pre(paste("Error:", conditionMessage(e))),
        footer = modalButton("Close")
      ))
    })
  })
  
  output$records_table <- renderDT({
    datatable(records_reactive(), options = list(pageLength = 5, scrollX = TRUE))
  })
  
  output$download_csv <- downloadHandler(
    filename = function() paste0("goldeneye_export_", Sys.Date(), ".csv"),
    content = function(file) {
      write.csv(records_reactive(), file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)