# this is stuff that needs to be run first and in order for S4 crap to work
methods::setOldClass("telemetry")
new.telemetry <- methods::setClass("telemetry", representation(info="list"), contains="data.frame")

methods::setOldClass("UD")
new.UD <- methods::setClass("UD", representation(info="list"), contains="list")

methods::setOldClass("ctmm")
new.ctmm <- methods::setClass("ctmm", representation(info="list"), contains="list")
