setClass("QcConfig",
         representation = list(
             file_paths = "character",
             groups = "numeric",
             group_names = "character",
             group_colors = "character"
         ))

QcConfig = function(file_paths,
                    groups = NULL,
                    group_names = NULL,
                    group_colors = NULL){
    if(is.null(groups)){
        groups = seq_along(file_paths)
    }
    if(is.null(group_names)){
        group_names = LETTERS[seq_along(unique(groups))]
    }
    if(is.null(group_colors)){
        group_colors = seqsetvis::safeBrew(length(group_names))
    }

    new("QcConfig",
        file_paths = file_paths,
        groups = groups,
        group_names = group_names,
        group_colors = group_colors
    )
}

setClass("QcConfigFeatures", contains = "QcConfig",
         representation = list(
             feature_load_FUN = "function",
             n_peaks = "numeric",
             min_fraction = "numeric",
             min_number = "numeric"
         ))

QcConfigFeatures = function(file_paths,
                            groups = NULL,
                            group_names = NULL,
                            group_colors = NULL,
                            n_peaks = 5e3,
                            min_fraction = 0,
                            min_number = 1){
    if(is.null(groups)){
        groups = seq_along(file_paths)
    }
    if(is.null(group_names)){
        group_names = LETTERS[seq_along(unique(groups))]
    }
    if(is.null(group_colors)){
        group_colors = seqsetvis::safeBrew(length(group_names))
    }

    new("QcConfigFeatures",
        file_paths = file_paths,
        groups = groups,
        group_names = group_names,
        group_colors = group_colors,
        n_peaks = n_peaks,
        min_fraction = min_fraction,
        min_number = min_number
    )
}

setClass("QcConfigSignal", contains = "QcConfig",
         representation = list(
             view_size = "numeric"
         ))

QcConfigSignal = function(file_paths,
                          groups = NULL,
                          group_names = NULL,
                          group_colors = NULL,
                          view_size = 3e3){
    if(is.null(groups)){
        groups = seq_along(file_paths)
    }
    if(is.null(group_names)){
        group_names = LETTERS[seq_along(unique(groups))]
    }
    if(is.null(group_colors)){
        group_colors = seqsetvis::safeBrew(length(group_names))
    }

    new("QcConfigSignal",
        file_paths = file_paths,
        groups = groups,
        group_names = group_names,
        group_colors = group_colors,
        view_size = view_size
    )
}

ob1 = QcConfig(LETTERS[1:5], groups = c(rep(1, 3), rep(2, 2)))
ob2 = QcConfigFeatures(LETTERS[1:5], groups = c(rep(1, 3), rep(2, 2)))
ob3 = QcConfigSignal(LETTERS[1:5], groups = c(rep(1, 3), rep(2, 2)))


setGeneric("QcColorMapping", function(object){standardGeneric("QcColorMapping")})
setMethod("QcColorMapping", c("QcConfig"), function(object){
    cols = object@group_colors
    names(cols) = as.character(object@group_names)
    cols
})


setGeneric("QcScaleColor", function(object){standardGeneric("QcScaleColor")})
setMethod("QcScaleColor", c("QcConfig"), function(object){
    cols = QcColorMapping(object)
    ggplot2::scale_color_manual(values = cols)
})

setGeneric("QcScaleFill", function(object){standardGeneric("QcScaleFill")})
setMethod("QcScaleFill", c("QcConfig"), function(object){
    cols = object@group_colors
    names(cols) = as.character(object@group_names)
    ggplot2::scale_fill_manual(values = cols)
})

QcColorMapping(ob3)
QcScaleFill(ob3)
