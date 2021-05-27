sys.source(paste("create_response_", response_model, ".R", sep = ""), envir = environment())

#create design matrix
Xy = data.frame(cbind(x_s, indic_T, y))
colnames(Xy) = c(paste0("x", 1 : p), "indic_T", "y")

#pull out yT, yC
yTs = Xy[Xy$indic_T == 1, "y"]
yCs = Xy[Xy$indic_T == 0, "y"]
