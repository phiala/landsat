slopeasp <- function (x, EWres, NSres, EWkernel, NSkernel, smoothing = 1, flatna=TRUE) 
{
    if(is(x, "SpatialGridDataFrame")) {
    	xmat <- t(as.matrix(x))
    }
    else {
       xmat <- as.matrix(x)
    }
    if (missing(EWres)) {
        if (is(x, "SpatialGridDataFrame")) {
            EWres <- slot(getGridTopology(dem), "cellsize")[1]
        }
        else {
            stop("EWres must be specified if x is not a SpatialGridDataFrame.\n")
        }
    }
    if (missing(NSres)) {
        if (is(x, "SpatialGridDataFrame")) {
            EWres <- slot(getGridTopology(dem), "cellsize")[2]
        }
        else {
            stop("NSres must be specified if x is not a SpatialGridDataFrame.\n")
        }
    }
    if (missing(EWkernel)) {
        EWkernel <- matrix(c(-1/8, 0, 1/8, -2/8, 0, 2/8, -1/8, 
            0, 1/8), ncol = 3, nrow = 3, byrow = TRUE)
    }
    EW.mat <- movingwindow(xmat, EWkernel)/EWres
    if (missing(NSkernel)) {
        NSkernel <- matrix(c(1/8, 2/8, 1/8, 0, 0, 0, -1/8, -2/8, 
            -1/8), ncol = 3, nrow = 3, byrow = TRUE)
    }
    NS.mat <- movingwindow(xmat, NSkernel)/NSres
    slope <- atan(sqrt(EW.mat^2 + NS.mat^2)/smoothing)
    slope <- (180/pi) * slope
    aspect <- 180 - (180/pi) * atan(NS.mat/EW.mat) + 90 * (EW.mat/abs(EW.mat))
    if(flatna) {
        aspect[slope == 0] <- NA
    } else {
        aspect[slope == 0] <- 0
    }

    akernel <- matrix(c(1,1,1,1,0,1,1,1,1), ncol = 3, nrow = 3, byrow = TRUE)
    smoothasp <- movingwindow(aspect, akernel)
    smoothcount <- movingwindow(ifelse(is.na(aspect), NA, 1), akernel)
    while(any(smoothcount) == 0) {
        smoothasp[smoothcount == 0] <- NA
        smoothcount[smoothcount == 0] <- NA
        smoothasp <- smoothasp/smoothcount
        smoothasp <- movingwindow(smoothasp, akernel)
        smoothcount <- movingwindow(ifelse(is.na(smoothasp), NA, 1), akernel)
    }
    smoothasp <- smoothasp/smoothcount
    aspect[slope == 0 & !is.na(slope)] <- smoothasp[slope == 0 & !is.na(slope)]


    if (is(x, "SpatialGridDataFrame")) {
        temp <- x
        temp@data[, 1] <- as.vector(t(aspect))
        aspect <- temp
        temp@data[, 1] <- as.vector(t(slope))
        slope <- temp
    }
    list(slope = slope, aspect = aspect)
}

