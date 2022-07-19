# ensure the solver is loaded using the ROI plugin mechanism; if
# the plugin is not available then offer to install it
check_ROI_solver <- function(solver) {
    if (solver != "lpsolve") {
        roi_plugin_name <- paste0("ROI.plugin.", solver)
        pkgload::check_suggested(roi_plugin_name, path = pkgload::inst("detectseparation"), version = "*")
        requireNamespace(roi_plugin_name, quietly = TRUE)
    }
}
