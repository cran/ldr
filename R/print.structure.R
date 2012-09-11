print.structure <-
function (x,...) 
{
	cat("\nLikelihood Ratio Test:\n");

	cat("Testing the hypotheses \n NH: Structure=", x$IC[1,1], "\n AH: Structure=", x$IC[2,1], "\n\n");

	print(x$LRT);

	cat("\nInformation Criteria:\n");

	print(x$IC);

	cat("----------\nNote on Structure: \niso -> isotropic \naniso -> anisotropic \nunstr -> unstructured\n")

	invisible(x)
}
