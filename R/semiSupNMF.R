
semiSupNMF <- function(
	A, k = 1L, 	init = NULL, mask = NULL,  check.k = TRUE,
	max.iter = 2000L, rel.tol = 1e-4, n.threads = 1L, trace = 1,
	verbose = 1L, show.warning = TRUE, inner.rel.tol = 1e-9
	) {

  	NNLM:::check.matrix(A, input.name = 'A');
	if (!is.double(A))
		storage.mode(A) <- 'double';
	n <- nrow(A);
	m <- ncol(A);

	init.mask <- NNLM:::reformat.input(init, mask, n, m, k);
	k <- init.mask$K;





	min.k <- min(dim(A));
	A.isNA <- is.na(A);
	A.anyNA <- any(A.isNA); # anyNA is depreciated in new version of R
	if (A.anyNA) {
	  stop("No NA allowed")
		# min.k <- min(min.k, ncol(A) - rowSums(A.isNA), nrow(A) - colSums(A.isNA));
		}
	rm(A.isNA);
	if (check.k && k > min.k )
		stop(paste("k larger than", min.k, "is not recommended, unless properly masked or regularized.
				Set check.k = FALSE if you want to skip this checking."));

	if (n.threads < 0L) n.threads <- 0L; # let openMP decide

	if (is.logical(verbose)) {
		verbose <- as.integer(verbose);
		}
	if (trace <= 0) {
		trace <- 999999L; # only compute error of the 1st and last iteration
		}

	run.time <- system.time(
		out <- c_semisupnnmf(A, as.integer(k),
			init.mask$Wi, init.mask$Hi, init.mask$Wm, init.mask$Hm,
			 as.integer(max.iter), as.double(rel.tol),
			as.integer(n.threads), as.integer(verbose), as.logical(show.warning),
			as.integer(1), as.double(inner.rel.tol),
			as.integer(trace))
		);
	names(out) <- c('W', 'H', 'mse', 'mkl', 'target.loss', 'average.epochs', 'n.iteration');
	out$mse <- as.vector(out$mse);
	out$mkl <- as.vector(out$mkl);
	out$target.loss <- as.vector(out$target.loss);
	out$average.epochs <- as.vector(out$average.epochs);

	# add row/col names back
	colnames(out$W) <- colnames(init.mask$Wi);
	rownames(out$H) <- rownames(init.mask$Hi);
	if (!is.null(rownames(A))) rownames(out$W) <- rownames(A);
	if (!is.null(colnames(A))) colnames(out$H) <- colnames(A);
	rm(init.mask);



	out$run.time <- run.time;
	out$options <- list(
		method = "semiSupervised_SCD",
		loss = "mkl",
		init = init,
		mask = mask,
		n.threads = n.threads,
		trace = trace,
		verbose = verbose,
		max.iter = max.iter,
		rel.tol = rel.tol,
		inner.max.iter = NA,
		inner.rel.tol = inner.rel.tol
		);
	out$call <- match.call();
	return(out);
	}
