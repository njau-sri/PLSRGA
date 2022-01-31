#' Gentype Construction

#' @param PCode: Marker genotype of parents original from diallel.
#' @param classifier: a label used for differentiated experimental design
#' @export
P2F1.diallel <- function(PCode, classifier = 'DF1'){

  # Purpose: This module used for constructing F1 genotype based on parent genotype of diallel
  #     ********Parameters**********
  # PCode: Marker genotype of parents original from diallel

  #     ********Format**************
  # First row is index of each parent
  # First column is name of each marker

  # Author: Wang Jinshe
  # E_mail: wangjinshe@gmail.com
  # Institute: NCSI (National Center for Soybean Improvement)
  # date: 2016/2/27

	nparent <- ncol(PCode)
	nmrk <- nrow(PCode)

	F1Code <- data.frame(ID = 1 : nmrk)
	cnms <- c()
	for (i in 1:nparent){
	 for (j in i:nparent){
		tmp <- matrix(paste(PCode[,i], PCode[, j], sep =""), ncol =1)
		F1Code <- cbind(F1Code, tmp)

		cnms <- c(cnms, paste('P', i, '.', j, sep = ''))
	 }
	}
	F1Code <- F1Code[, -1]

	colnames(F1Code) <- cnms
	if (length(classifier) == 1)
		Classifier <- matrix(c(rep(classifier, ncol(F1Code))), nrow = 1)
	else if (length(classifier) > 1 && length(classifier) < ncol(F1Code))
		stop('Length of classifier is different with the number of F1 crosses, please check it!')
	else
		Classifier <- classifier

	if (length(Classifier) == ncol(F1Code)){
		colnames(Classifier) <- cnms
		F1Code <- rbind(Classifier, F1Code)
		rownames(F1Code) <- c('class', rownames(PCode))
	}

	return(F1Code)
}

#' Gentype Construction

#' @param PCode: Marker genotype of parents original from NCII
#' @param classifier: a label used for differentiated experimental design
#' @export
P2F1.NCII <- function(PCode, classifier = 'NCF1'){

  # Purpose: This module used for constructing F1 genotype based on parent genotype of diallel
  #     ********Parameters**********
  # PCode: Marker genotype of parents original from diallel

  #     ********Format**************
  # First row is index of each parent
  # First column is name of each marker

  # Author: Wang Jinshe
  # E_mail: wangjinshe@gmail.com
  # Institute: NCSI (National Center for Soybean Improvement)
  # date: 2016/2/27

	group <- PCode[1, ]
	PCode <- PCode[-1, ]
	ID <- unlist(unique(c(group)))
	nmale <- sum(group == ID[1])
	nfemale <- sum(group == ID[2])

	NCF1 <- NULL
	cnms <- c()
	for (i in 1:nmale){
		for (j in (nmale + 1):(nmale + nfemale)){
			tmp <- matrix(paste(PCode[, i], PCode[, j], sep = ''))
			NCF1 <- cbind(NCF1, tmp)
			cnms <- c(cnms, paste('P', i, '.', j, sep = ''))
		}
	}

	for (i in 1:(nmale + nfemale)){
		tmp <- matrix(paste(PCode[, i], PCode[, i], sep = ''))
		NCF1 <- cbind(NCF1, tmp)
		cnms <- c(cnms, paste('P', i, '.', i, sep = ''))
	}
	colnames(NCF1) <- cnms
	if (length(classifier) == 1)
		Classifier <- matrix(c(rep(classifier, ncol(NCF1))), nrow = 1)
	else if (length(classifier) > 1 && length(classifier) < ncol(F1Code))
		stop('Length of classifier is different with the number of F1 crosses, please check it!')
	else
		Classifier <- classifier
	#print(Classifier)
	if (length(Classifier) == ncol(NCF1)){
		colnames(Classifier) <- cnms
		NCF1 <- rbind(Classifier, NCF1)
		#print(nrow(NCF1))
		tmp <- c('class', rownames(PCode)[-1])
		#print(length(tmp))
		rownames(NCF1) <- c('class', rownames(PCode))
	}

	return(NCF1)
}

#' Fitness function

#' @param chr: GA Chromosome
#' @param parent: The DataSet
#' @param tr: Train dataset
#' @param te: test dataset
#' @param res: residuals
reg.fitness <- function(chr, parent, tr, te, res){

  # Purpose: fitness function of partial least squares (pls)
  #**************--------parameters-------**************
  # chr: chromosome
  # parent: the dataset
  # tr: train dataset
  # te: test dataset
  # res: residuals

  # Author: Wang Jinshe
  # E_mail: wangjinshe@gmail.com
  # Institute: NCSI (National Center for Soybean Improvement)
  # date: 2016/2/27

  try <- parent$data$dependent
  trd <- data.frame(parent$data$data[, as.numeric(chr)])
  trm <- lm(try ~ ., data = trd)
  tey <- parent$data$dependent[te]
  ted <- data.frame(parent$data$data[te, as.numeric(chr)])
  cor(predict(trm), try)^2
}

reg.fitness.pls <- function(chr, parent, tr, te, res){

  try <- parent$data$dependent
  tmp0  <- mrk[as.numeric(chr), ]
  tmp0 <- apply(tmp0, 1, F1Construct)
  mrkf1 <- do.call(rbind, tmp0)
  trd <- as.data.frame(t(mrkf1))
  pls.out <- plsr(try ~ ., method = 'simpls', validation = 'CV', scale = F, data=trd)
  expvar <- as.matrix(explvar(pls.out))
  nn <- length(expvar)
  dif <- expvar[1 : (nn - 1)] - expvar[2 : nn]
  kk <- which(dif < 0)
  if (length(kk) == 0){
    best <- nn
  }
  else {
    best <- kk[1]
  }
  tey <- parent$data$dependent[te]
  ted <- as.data.frame(t(mrkf1))[te, ]
  cor(predict(pls.out, ncomp=best, newdata=trd), try)^2
}

#' Single Marker Analysis

#' @param phe: Marker genotype of parents original from diallel or NCII design
#' @param mrk: column: Individuals, row: Markers
#' @export
single.MrkAna <- function(phen, mrk){
 # Purpose: single marker analysis
 # mrk: column: Individuals
 #      row: Markers

	nmrk <- nrow(mrk)
	pvalue <- matrix(0, nrow = nmrk, ncol =1)
	for (i in 1:nmrk){
		lm.out <- lm(phen ~ unlist(mrk[i, ]))
		if (i %% 500 == 0) print('+')
		fstat <- summary(lm.out)$fstatistic
		pp <- 1 - pf(fstat[1], fstat[2], fstat[3])
		pvalue[i] <- pp
	}
	return(pvalue)
}

#' F1 genotype Construction

#' @param mrk: genotypes of each parent for diallel design, First row including parent's ID, First column including marker's ID.
F1Construct <- function(mrk){

	# Purpose: Construction of design matrix for F1
	#
	#***********---Parameters---**************
	# mrk: genotypes of each parent for diallel design
	#      First row including parent's ID
	#      First column including marker's ID

	# Author: Wang Jinshe
	# Email: wangjinshe@gmail.com
	# Institute: NCSI (National Center for Soybean Improvement)
	# Date: 2016/02/29

	allel <- matrix(unique(c(mrk)), nrow = 1)
	mrkname <- rownames(mrk)
	mrk <- t(as.matrix(mrk))
    #nmrk <- nrow(mrk)
	nP <- length(mrk)
	nF1 <- nP * (nP + 1) / 2

	nallel <- length(allel)
	nn <- nallel * (nallel + 1) / 2
	F1Code <- matrix(0, nn, nF1)
	flag = 1
	allel_name <- c()
	for (i in 1:nP){
		for (j in i:nP){
			id1 <- which(allel == mrk[i])
			id2 <- which(allel == mrk[j])
			if (id1 == id2){
				F1Code[id1, flag] <- 1
				allel_name <- c(allel_name, paste('A', '.',
					allel[id1], '.', allel[id1], sep = ''))
			} else
			{
				F1Code[id1, flag] <- 0.5
				F1Code[id2, flag] <- 0.5
				if (id1 < id2){
					idd <- id1 / 2 * (2 * nallel - id1 + 1) + id2 - id1
				}
				else{
					idd <- id2 / 2 * (2 * nallel - id2 + 1) + id1 -id2
				}
				F1Code[idd, flag] <- 1
				allel_name <- c('D', paste(mrkname, '.',
					allel[id1], '.', allel[id2], sep = ''))
			}
			flag = flag + 1
		}
	}
	add_names <- paste('A', allel, sep = "")
	dom_names <- ""
	if (nallel >2){
		for (i in 1:(nallel-1)){
			for (j in (i + 1):nallel){
				dom_names <- c(dom_names, paste('D', allel[i], '.', allel[j], sep = ''))
			}
		}
		dom_names <- dom_names[-1]
	}
	else if (nallel == 2){
		dom_names <- paste('D', allel[1], '.', allel[2], sep = '')
	}
	ad_names <- c(add_names, dom_names)
	rownames(F1Code) <- ad_names
    F1Code <- as.data.frame(F1Code)
	return(F1Code)
}

#' partial least square regression

#' @param phen: the phenotype of hybrid trait, 1 column
#' @param mrk: the design matrix, output of module F1Construct
#' export
GA.PLS <- function(phen, mrk){
 # Purpose: partial least square regression (PLSR)
 # phen: The phenotype
 # mrk: the design matrix, output of module F1Construct

 # Author: Wang Jinshe
 # Email: wangjinshe@gmail.com
 # Institute: NCSI (National Center for Soybean Improvement)
 # Date: 2016/03/06

 pls.out <- plsr(phen ~ mrk, method = 'simpls', ncomp = length(phen)-3, validation = 'LOO', scale = F)
 expvar <- as.matrix(explvar(pls.out))
 nn <- length(expvar)
 dif <- expvar[1 : (nn - 1)] - expvar[2 : nn]
 kk <- which(dif < 0)
 if (length(kk) == 0){
	best <- nn
 }
 else {
	best <- kk[1]
 }
 g.effect <- coefficients(pls.out, ncomp = best, intercept = TRUE)
 R.square <- (cor(as.matrix(predict(pls.out, ncomp = best)), phen))^2
 list(geffect = g.effect, Rsquare = R.square)
}

#' Estimation of R square and genetic variance of each loci

#' @param phe: the phenotype of hybrid trait, 1 column
#' @param mrk: the design matrix, output of module F1Construct
#' @param geffect: additive and dominance effect matrix
#' export
R.Square <- function(mrk, geffect, Phen, Rsquare){
 # Purpose: Estimation genetic variance of each locus
 # object: output of plsr module

 # Author: Wang Jinshe
 # Email: wangjinshe@gmail.com
 # Institute: NCSI
 # Date: 2016/03/07

 nmrk <- length(mrk)
 id <- 1
 intercept <- geffect[1]
 geffect <- geffect[2:length(geffect)]
 r2 <- NULL
 adv <- matrix(0, nrow = nmrk, ncol = 2)
 for (i in 1:nmrk){
	nn <- nrow(mrk[[i]])
	tmp <- matrix(geffect[id:(id + nn -1)], nrow = 1)
	tmp1 <- mrk[[i]][substr(rownames(mrk[[i]]), 1, 1) == 'A', ]
	tmp2 <- mrk[[i]][substr(rownames(mrk[[i]]), 1, 1) == 'D', ]

	id1 <- nrow(tmp1)
	cc <- NULL
	for (j in 1:id1){
		cc <- c(cc, as.numeric(tmp1[j, ] * tmp[1:id1][j] ))
	}
	vA <- var(cc)
	cc <- NULL
	for (j in 1:nrow(tmp2)){
		cc <- c(cc, as.numeric(tmp2[j, ] * tmp[(id1 + 1) : nn]))
	}
	vD <- var(cc)
	adv[i, 1] <- vA
	adv[i, 2] <- vD
	pred <- matrix(geffect[id:(id + nn -1)], nrow = 1) %*% as.matrix(mrk[[i]]) + intercept
	r2 <- c(r2, cor(matrix(c(pred), ncol = 1), matrix(Phen))^2)

	id <- id + nn
 }
 percent <- r2 / sum(r2)
 R2 <- matrix(as.vector(percent) * as.vector(Rsquare), ncol = 1)
 rownames(adv) <- names(mrk)
 adR2 <- cbind(adv, R2)
 colnames(adR2) <- c('V_A', 'V_D', 'R2')
 return(adR2)
}

#' Estimation of degree of dominance of each loci

#' @param addmatrix: addtive matrix
#' @param dommatrix: dominance matrix
#' export
DD <- function(){
    addmatrix <- read.csv(file.choose())
    dommatrix <- read.csv(file.choose())

    tmp <- do.call(rbind, strsplit(as.character(addmatrix[,1]), split = "@"))[, 1]
    selmrk <- unique(tmp)
    tmp1 <- do.call(rbind, strsplit(as.character(dommatrix[, 1]), split = "@"))[, 1]
    aveADEst <- function(sel){
        subadd <- addmatrix[tmp == sel, ]
        prob_add <- rowSums(subadd[, 3:ncol(subadd)])/sum(subadd[, 3:ncol(subadd)])
        ave_add <- matrix(as.numeric(subadd[, 2]), nrow = 1) %*% matrix(prob_add, ncol = 1)
        subdom <- dommatrix[tmp1 == sel, ]
        prob_dom <- rowSums(subdom[, 3:ncol(subdom)])/sum(subdom[, 3:ncol(subdom)])
        ave_dom <- matrix(as.numeric(subdom[, 2]), nrow = 1) %*% matrix(prob_dom, ncol = 1)
        return(c(ave_add, ave_dom))
    }
    ave_ad <- lapply(selmrk, aveADEst)
}

#' Estimation genetic effect of each loci

#' @param add1: additive effect1
#' @param add2: additive effect2
#' @param dom: dominance
#' export
geffect <- function(add1, add2, dom){
	#--------Parameters------------
	# add1, add2: additive effect
	# dom: dominance effect

	est <- 0.5 * add1 + 0.5 * add2 + dom
	return(est)
}

#' Identify duplicated molecular marker

#' @param mrk: The original parent's genotypes
#' @param selmrk: output of module Sel.mrk
#' export
Dup.find <- function(mrk, selmrk){

  # Purpose: Identify duplicated molecular marker
  # Augument:
  #    mrk: The original parent's genotypes
  #    selmrk: output of module Sel.mrk

  # Author: Wang Jinshe
  # Email: wangjinshe@gmail.com
  # Institute: NCSI (National Center for Soybean Improvement)
  # Date: 2016/03/28

	mrk <- as.matrix(mrk)
	nsel <- length(selmrk)
	selmrkcode <- mrk[rownames(mrk) %in% selmrk, ]
	dup.pos <- list()
    for (i in 1:nsel){
        cat(i, ' th marker\n', sep = "")
		id <- which(apply(mrk, 1, function(x) identical(x, selmrkcode[i, 1:ncol(mrk)])))
		dup.pos[[i]] <- mrk[id, ]

	}
	names(dup.pos) <- selmrk
	return(dup.pos)
}

#' The module is used for constructiong allele matrix

#' @param mrk: marker genotype of each parent
#' @param selmrk: F1 design matrix of identification markers
#' @param est.para: a list output from est.para module, which include genetic effect and R-square estimation
#' @param design: a string, 'NCII' or 'diallel'.
#' export
Allele.matrix <- function(mrk, selmrk, est.para, design = "diallel"){

	# Purpose: The module is used for constructiong allele matrix
	# mrk: marker genotype of each parent
	# selmrk: F1 design matrix of identification markers
	# est.para: a list output from est.para module, which include genetic effect and R-square estimation

	# Author: Wang Jinshe
	# Email: wangjinshe@gmail.com
	# Institute: NCSI (National Center for Soybean Improvement)
	# Date: 2016/03/18

	np <- ncol(mrk)
	#nmrk <- length(selmrk)
	if (design == "NCII"){
		Mean <- est.para$geffect[1]
		geffect <- est.para$geffect[-1]
		id <- mrk[1, ]
		nr <- sum(id == 1)
		nc <- sum(id == 2)
		mrk <- mrk[rownames(mrk) %in% selmrk, ]
		mrk <- rbind(id, mrk)

		f1.mrk <- F1Construct(mrk)
		nn <- sapply(f1.mrk, nrow)
		mrk1 <- do.call(rbind, f1.mrk)
		mrkname <- paste(rep(selmrk, nn), rownames(mrk1), sep ="@")
		rownames(mrk1) <- mrkname

		add.matrix <- mrk1[, (ncol(mrk1) - np +1):ncol(mrk1)]
		id <- rowSums(add.matrix) != 0
		add.matrix <- add.matrix[id, ]
		add.matrix <- cbind(matrix(geffect[id], ncol = 1), add.matrix)

		tmp <- mrk1[, 1:(ncol(mrk1) - np)]
		id <- apply(tmp, 1, function(x) !any(x == 0.5))
		tmp <- tmp[id, ]
		dom <- geffect[id][rowSums(tmp) != 0]
		tmp <- tmp[rowSums(tmp) != 0, ]
		dom.matrix <- tmp
		dom.matrix <- cbind(matrix(dom, ncol = 1), dom.matrix)
		tmp <- expand.grid(paste('C', 1 : nc, sep = ""), paste("R", 1 : nr, sep = ""))
		tmp <- paste(tmp[, 1], tmp[, 2], sep = ".")
		colnames(dom.matrix) <- c("Estimation", tmp)
	}
	else if (design == "diallel"){
	    mrk <- mrk[rownames(mrk) %in%selmrk, ]
		f1.mrk <- apply(mrk, 1, F1Construct)
		nn <- sapply(mrk, nrow)
		mrk1 <- do.call(rbid, f1.mrk)
		mrkname <- paste(rep(selmrk, nn), rownames(mrk1), sep = "@")
		rownames(mrk1) <- mrkname
		id <- rowSums(mrk1) != 0
		add.matrix <- mrk1[, (ncol(mrk1) - np + 1):ncol(mrk1)]

		tmp <- mrk1[, 1:(ncol(mrk1) - np)]
		id <- apply(tmp, 1, function(x) !any(x == 0.5))
		tmp <- tmp[id, ]
		tmp <- tmp[rowSums(tmp) != 0, ]
		dom.matrix <- tmp
		tmp <- expand.grid(paste("P", 1:np, sep = ""), paste("P", 1:np, sep = ""))
		tmp <- paste(tmp[, 2], tmp[, 1], sep = ".")
		colnames(dom.matrix) <- tmp
	}
	list(addMatrix = add.matrix, domMatrix = dom.matrix)
}

#' The module is used for additive, dominance and R-square estimation

#' @param phe: phenotype
#' @param mrk: marker genotype of each parents
#' @param ga.out: BigBang object original from configBB.VarSelMisc (package: galgo)
#' export
EST.para <- function(phe, mrk, ga.out){
	# Purpose: additive, dominance and R-square estimation
	# phe: phenotype
	# mrk: marker genotype of each parents
	# ga.out: BigBang object original from configBB.VarSelMisc (package: galgo)

	# Author: Wang jinshe
	# Email: wangjinshe@gmail.com
	# Institute: NCSI (National Center for Soybean Improvement)
	# Date: date()
	if (!require(pls))
	    library(pls)
	if (sum(class(ga.out) == 'BigBang') != 1)
	    stop('The third argument for this module should be a BigBang object')

	selmrk <- Sel.mrk(ga.out)
	mrkname <- rownames(mrk)
	mrk <- rbind(mrk[1, ], mrk[mrkname %in% selmrk, ])
	mrk <- F1Construct(mrk)
	nn <- sapply(mrk, nrow)
	mrk1 <- do.call(rbind, mrk)
	mrkname <- paste(rep(selmrk, nn), rownames(mrk1), sep = "@")
	rownames(mrk1) <- mrkname
	mrk1 <- t(mrk1)
	gapls.out <- GA.PLS(phe, mrk1)
	print('****Finished****')
	R2 <- R.Square(mrk, gapls.out$geffect, phe, gapls.out$Rsquare)

	list(geffect = gapls.out$geffect, varR2 = R2)
}

GEffect <- function(add1, add2, dom){
  # If interested locus is homozygous, then add1 = add2 = allele effect and dom = 0
  # If the locus is heterogous, then add1 = allele 1 effect, add2 = allele 2 effect and
  #dom = dominance effect

  gf = 2 * add1 + 2 * add2 + dom
  return(gf)
}

gpred <- function(add, id_add, dom, id_dom, index, aname){
  # Compute the prediction of present estimation

  mrk <- NULL
  allele <- NULL
  for (i in 1:length(aname)){
    tmp <- strsplit(as.character(aname[i]), split = "@")[[1]]
    mrk <- c(mrk, tmp[1])
    allele <- c(allele, tmp[2])
  }
  nF1 <- ncol(index)
  pred_gf <- NULL
  m_a_a <- NULL
  for (i in 1:nF1){
    if (any(id_dom[, i] != 0)){
      E_dom <- dom[as.logical(id_dom[, i])]
    }
    else {
      E_dom <- 0
    }
    id1 <- id_add[, index[1, i]]
    id2 <- id_add[, index[2, i]]
    allele1 <- add[as.logical(id1)]
    allele2 <- add[as.logical(id2)]
    tmp <- matrix(c(GEffect(allele1, allele2, E_dom),
                    GEffect(allele1, allele1, 0),
                    GEffect(allele2, allele2, 0)),
                  nrow = 1)
    pred_gf <- rbind(pred_gf, tmp)
    pred_homo1 <-
    tmp <- matrix(c(mrk[as.logical(id1)],
                    paste("P", index[1, i], sep = ""),
                    allele[as.logical(id1)],
                    paste("P", index[2, i], sep = ""),
                    allele[as.logical(id2)]
                    ),
                  nrow = 1)
    m_a_a <- rbind(m_a_a, tmp)
  }
  mrk_vale = data.frame(mrk = m_a_a[, 1],
                        P1 = m_a_a[, 2],
                        A1 = m_a_a[, 3],
                        P2 = m_a_a[, 4],
                        A2 = m_a_a[, 5],
                        Heter = pred_gf[, 1],
                        Homo1 = pred_gf[, 2],
                        Homo2 = pred_gf[, 3]
                        )

  return(mrk_vale)
}

filt <- function(MrkVale){
  Dom_value <- unique(MrkVale$Heter)
  A1 <- NULL
  A2 <- NULL
  for (i in 1:length(Dom_value)){
    tmp <- MrkVale[which(MrkVale$Heter == Dom_value[i]), ]
    A1 <- c(A1, as.character(tmp$A1)[1])
    A2 <- c(A2, as.character(tmp$A2)[1])
  }
  AD_value <- data.frame(Value = Dom_value, A1 = A1, A2 = A2)
  return(AD_value)
}

Improv <- function(MrkValue, AD_value){
  conn <- file(file.choose(), 'a')
  cat(paste('Locus ',
            MrkValue[1,]$mrk,
            ' Alternative Improvement Programs are listed as follow: \n',
            sep=''), file = conn)
  cat('**********************************************\n', file = conn)
  nF1 <- nrow(MrkValue)
  for (i in 1:nF1){
    heter_vale <- MrkValue[i, 6]
    A1 <- MrkValue[i, 3]
    A2 <- MrkValue[i, 5]
    tmp <- AD_value[AD_value$Value > heter_vale, ]
    if (length(tmp$Value) > 1){
      cat('------------------------\n', file = conn)
      cat(paste(MrkValue[i, ]$P1, '/', MrkValue[i, ]$P2, '\n', sep = ''), file = conn)
      cat(paste(MrkValue[i, ]$A1, '/',
                  MrkValue[i, ]$A2,
                  '(',
                  MrkValue[i, ]$Heter,
                  ')', ' -->> \n',
                  sep = ''), file = conn)
      for (j in 1:length(tmp$A1)){
        cat(paste(tmp[j, ]$A1,
                  '/',
                  tmp[j, ]$A2,
                  '  Value = ',
                  tmp[j, ]$Value,
                  '\n',
                  sep = ''), file = conn)
      }
      cat('------------------------\n', file = conn)
    }
  }
  cat('**********************************************\n', file = conn)
  cat('***************--The End--*********************\n', file = conn)
  close(conn)
}

MainImprov <- function(add, dom){
  index <- dom[1:2, 3:ncol(dom)]
  dom <- dom[-(1:2), ]
  amrk <- NULL
  for(i in 1:nrow(add)){
    amrk <- c(amrk, strsplit(as.character(add[i, 1]), split = '@')[[1]][1])
  }
  dmrk <- NULL
  for (j in 1:nrow(dom)){
    dmrk <- c(dmrk, strsplit(as.character(dom[j, 1]), split = '@')[[1]][1])
  }
  mrk <- unique(amrk)
  for (i in 1:length(mrk)){
    cat(paste(mrk[i],'\n'))
    id <- amrk == mrk[i]
    a_effect <- add[id, 2]
    a_index <- add[id, 3:ncol(add)]
    a_name <- add[id, 1]
    id <- dmrk == mrk[i]
    d_effect <- dom[id, 2]
    d_index <- dom[id, 3:ncol(dom)]
    geffect_pred <- gpred(a_effect, a_index, d_effect, d_index, index, a_name)
    Filt <- filt(geffect_pred)
    print(geffect_pred)
    print(Filt)
    Improve <- Improv(geffect_pred, Filt)
  }
}

ImpInfoExact <- function(...){

    # Input the output of module MainImprov()
    #    a .txt file should be selected

    #imp <- readLines(file.choose())

    path <- file.choose()
    ddd <- basename(path)
    if (strsplit(as.character(ddd), split = ".", fixed = T)[[1]][2] != "txt")
        stop("The extension of input file should be '.txt', please check it!")

    imp <- readLines(path)
    ID <- grep("^Locus", imp)
    nmrk <- length(ID)
    id <- cbind(ID[1:(nmrk - 1)], ID[2:nmrk] - 1)
    id <- rbind(id, matrix(c(ID[nmrk], length(imp)), nrow = 1))
    id[, 1] <- id[, 1] + 2 # Index for each marker

    # Exact Locus
    locus <- sapply(imp[ID], function(x) strsplit(x, split = " ")[[1]][2])
    names(locus) <- locus

    # Exact cross
    #cross_id <- grep("^P[0-9]/P[0-9]", imp, perl = T)
    #cross <- unique(imp[cross_id])
    idd <- data.frame(Marker = locus, Start = id[, 1], End = id[, 2])

    impreshape <- function(x){
        #print(as.numeric(x[3]))
        tmp <- imp[as.numeric(x[2]):(as.numeric(x[3]) - 2)]
        lowerline <- grep("------------------------", tmp, perl = T)
        idcross <- cbind(lowerline[seq(1, length(lowerline) -1, 2)], lowerline[seq(1, length(lowerline) - 1, 2) + 1])
        cross_mrk <- tmp[idcross[, 1] + 1]

        PresentGValue <- tmp[idcross[, 1] + 2]
        PresentAllele <- do.call(rbind, strsplit(as.character(PresentGValue), split = "\\("))[, 1]
        id3 <- matrix(unlist(gregexpr("[()]", PresentGValue, perl = T)), ncol = 2, byrow = T)
        PresentGValue <- substr(PresentGValue, id3[, 1] + 1, id3[, 2] - 1)
        PresentGValue <- as.numeric(PresentGValue)

        PredGValue <- function(x){
            pred <- tmp[(x[1] + 3):(x[2] - 1)]
            tt <- do.call(rbind, strsplit(pred, split = " "))
            pred <- data.frame(Allele = tt[, 1], Value = tt[, 5])
        }

        pred <- apply(idcross, 1, PredGValue)
        nn <- sapply(pred, nrow)
        cross_mrk <- rep(cross_mrk, nn)
        PresentGValue <- rep(PresentGValue, nn)
        pred <- do.call(rbind, pred)
        PresentAllele <- rep(PresentAllele, nn)
        #PresentGValue <- as.numeric(levels(PresentGValue))[as.integer(PresentGValue)]
        #print(PresentGValue)
        predvalue <- pred[, 2]
        predvalue <- as.numeric(levels(predvalue))[as.integer(predvalue)]
        Result <- data.frame(Cross = cross_mrk, PresAllele = PresentAllele, PredAllele = pred[, 1], PresGValue = PresentGValue, PredValue = predvalue)
        #Result$Marker <- rep(x[1], nrow(Result))
        return(Result)
    }

    ImpInfo <- apply(idd, 1, impreshape)
    ImpInfo <- do.call(rbind, ImpInfo)
    return(ImpInfo)
}

ImprovPlanReduc <- function(impinfo){

    cross <- unique(impinfo$Cross)
    ncross <- length(cross)
    Result <- NULL
    ss <- NULL
    for (i in cross){
        print(i)
        tmp <- impinfo[impinfo$Cross == i, ]
        mrks <- do.call(rbind, strsplit(rownames(tmp), split = ".", fixed = TRUE))[, 1]
        smrk <- unique(mrks)
        ttt <- NULL
        #print(class(tmp))
        for (j in smrk){
            #cat(paste(j, "...", sep = "  "))
            tmp1 <- tmp[mrks == j, ]
            #print(names(tmp1))
            idmin <- which(tmp1$PredValue == min(tmp1$PredValue))
            idmax <- which(tmp1$PredValue == max(tmp1$PredValue))
            min_max_pred <- rbind(tmp1[idmin, ], tmp1[idmax, ])
            #print(tmp1)
            #print(idmin)
            #print(idmax)
            #print(min_max_pred)
            ttt <- rbind(ttt, data.frame(Cross = min_max_pred$Cross[1], PresAllele = min_max_pred$PresAllele[1], minPredAllele = min_max_pred$PredAllele[1], maxPredAllele = min_max_pred$PredAllele[2], PresGValue = min_max_pred$PresGValue[1], minPredValue = min_max_pred$PredValue[1], maxPredValue = min_max_pred$PredValue[2]))
        }
        ss <- c(ss, unlist(smrk))

        ttt$CrossMinPred <- rep(sum(ttt$minPredValue), nrow(ttt))
        ttt$CrossMaxPred <- rep(sum(ttt$maxPredValue), nrow(ttt))
        Result <- rbind(Result, ttt)
    }
    Result$mrk <- ss
    return(Result)
}

#' The module used for Marker selection

#' @param BigBang: output of GAPLS
#' @export
Sel.mrk <- function(BigBang){

 # Purpose: Marker selection
 # BigBang: BigBang Object

 # Author: Wang Jinshe
 # Email: wangjinshe@gmail.com
 # Institute: NCSI (National Center for Soybean Improvement)
 # Date: 2016/03/01

	OK <- filterSolution(BigBang, "solutions", TRUE)
	x <- geneFrequency(BigBang, "solutions", OK, cutoff = -1, gene.names = FALSE)
	r <- unique(x)
	v <- sapply(r, function(i) sum(x >= i))
	rcol <- c(cut(1:50, breaks = 8, label = FALSE), 0)
	n.sel <- sum(!is.na(rcol[v]))

	x <- geneFrequency(BigBang, gene.names = TRUE)
	tmp <- sort(x, decreasing = TRUE)
	sel.mrk <- names(tmp[1:n.sel])
	sel.mrk <- gsub(pattern = " ", replacement = "", sel.mrk)
	sel.mrk <- as.matrix(sel.mrk)

	return(sel.mrk)
}

#' The module used for QTL idenfication and genetic parameter estimation

#' @param filepath: a string, The path to the folder where the original data and analysis results are stored
#' @param mrkfile: Parental genotype file in CSV format
#' @param traitfile: F1 phenotyp file in csv format
#' @export
PLSRGA <- function(filepath, mrkfile, traitfile,
		Exp.Des = "Diallel", singmrk = FALSE, para = TRUE,
		chromSize = 5, max_gnr = 1000, pop_size = 50, goal_fit = 0.95, save_freq = 50){

 # Purpose: partial least square regression (PLSR)
 # filepath: a string, The path to the folder where the original data and analysis results are stored
 # mrkfile: parental genotype file in csv format，Example: 'D:/data/marker.csv'
 # traitfile: F1 phenotyp file in csv format; Example: ‘D:/data/phenotype.csv’

 # Author: Wang Jinshe
 # Email: wangjinshe@gmail.com
 # Institute: NCSI (National Center for Soybean Improvement)
 # Date: 2016/03/06


# load data
  # -----------------++++ Structure of Marker genotype of parents ++++--------------------
    # The data source is a text file with tab delimited or csv file with commas separated
    # according your favoriate. The expected file format is marker genotypes in rows and
    # samples in columns. The first columns must be marker names identifier, accession number,
    # or anything to distinguish uniquely the markers. The first row must contain the sample
    # names (identifier of parents using for diallel or NCII), again unique values. In NCII design,
    # the second row is the class description for each parents group, in other words, the second
    # row contain the identifier that were used for distinguish male and female parents. Example:
    #			P1	P2	P3	P4	P5	P6	P7	P8
	#	Mrk1	3	3	3	3	3	1	1	3
	#	Mrk2	4	4	4	4	4	1	1	4
	#	Mrk3	1	1	1	1	1	4	4	1
	#	Mrk4	3	3	3	3	3	2	2	3
	#	Mrk5	3	3	3	3	3	4	4	3
	#	Mrk6	3	3	3	3	3	2	2	3

	# ... ...
  #	----------------------++++ Structure of phenotype ++++--------------------------------
    # Trait were ordered by column, in other words, one trait was included in one column. The
    # most important issue is the rank of each row.
    # In marker dataset, the column ranked as m1, m2, m3,... f1, f2, f3,..., 'mi' and 'fi' denote
    # male and female parents respectively. Phenotype should be ordered as:
    #        code	phe
	#		P1*P2	110.4055258
	#		P1*P3	382.5267308
	#		P1*P4	430.4869853
	#		P1*P5	393.5830995


#------------------------------- Install required packages --------------------------------
# Install required packages pls and galgo

	# rm(list = ls())
	if (!library(pls, logical.return=T)){
		chooseCRANmirror()
		install.packages('pls')
	}
	if (!library(galgo, logical.return=T)){
		chooseCRANmirror()
		install.packages(c("R.oo", "MASS", "class", "e1071", "rpart", "nnet", "randomForest"))
		install.packages('galgo')
	}
	library(pls)
	library(galgo)

	if (is.null(filepath)){
		stop("The file path must be scheduled!")
	}
	if (is.null(mrkfile)){
		stop("The marker gentoype file name must be scheduled!")
	}
	if (is.null(traitfile)){
		stop("The trait file name must be scheduled!")
	}

	setwd(filepath)
	mrk <- read.csv(mrkfile, row.names = 1, header=T)
	phe <- read.csv(traitfile, row.names = 1, header=T)
	assign("mrk", mrk, envir=.GlobalEnv)

# create directory which is used for storage the result file
	dirs <- dir()
	if (sum(dirs == "Result") == 0){
		dir.create("Result")
	}
	else {
		cat("The folder 'Result' already exists, please check it and keep empty!")
	}

# checking the data type, structure .et al

# Construction of F1 structure matrix
	if (toupper(Exp.Des) == 'DIALLEL'){
		tmp <- P2F1.diallel(mrk)
	}
	else if (toupper(Exp.Des) == 'NCII'){
		tmp <- P2F1.NCII(mrk)
	}
# Single marker analysis to filter the dataset
	if (singmrk){
		id <- tmp[1, ]
		tmp1 <- tmp[-1, ]
		pval <- single.MrkAna(phe, tmp1)
		id1 <- pval < 0.05
		tmp1 <- tmp1[id1, ]
		tmp <- rbind(id, tmp1)
	}
	tmp <- unique(tmp)

	label <- rownames(tmp)
	tmp <- cbind(label, tmp)
	write.table(tmp, 'tmp.txt', row.names=F, quote = FALSE, sep = "\t")

# Loci identification using GA based the R package galgo
	if (para){
		bb <- configBB.VarSelMisc(file = "tmp.txt", chromosomeSize = chromSize, niches = 1, populationSize = pop_size,
		maxGenerations = max_gnr, goalFitness = goal_fit, saveVariable = "bb", saveFrequency = save_freq,
		saveFile = "GAoutput.parallel.Rdata", fitnessFunc = reg.fitness.pls,
		callEnhancerFunc = function(chr, parentBB){robustGeneBackwardElimination(chr, parentBB, result = "shortest")})

		bb$data$dependent <- phe[, 1]
		assignParallelFile(bb)
		blast(bb)
	}
	else {
		bb <- configBB.VarSelMisc(file = "tmp.txt", chromosomeSize = chromSize, niches = 1, populationSize = pop_size,
		maxGenerations = max_gnr, goalFitness = goal_fit, saveVariable = "bb", saveFrequency = save_freq,
		saveFile = "GAoutput.Rdata", fitnessFunc = reg.fitness.pls)
		bb$data$dependent <- phe[, 1]
		blast(bb)
	}

# Estimation of genetic parameters using partial least squares (pls) based on the R package pls
	if (toupper(Exp.Des) == "DIALLEL"){
		selmrk <- Sel.mrk(bb)
		mrkname <- rownames(mrk)
		submrk <- mrk[mrkname %in% selmrk, ]
		subGenotype <- apply(submrk, 1, F1Construct)
		nn <- sapply(subGenotype, nrow)
		mrk1 <- do.call(rbind, subGenotype)
		mrk1 <- t(mrk1)
		gapls.out <- GA.PLS(phe[,1], mrk1)
		R2 <- R.Square(subGenotype, gapls.out$geffect, phe[,1], gapls.out$Rsquare)
		est.para <- list(geffect = gapls.out$geffect, varR2 = R2)
	}
	else if (toupper(Exp.Des) == "NCII"){
		est.para <- EST.para(phe, mrk, bb)
	}

# Collating the results and additive, dominance effect and R square of each locus were exported
	geffect <- est.para$geffect
	varR2 <- est.para$varR2
	write.csv(geffect, "Result/geffect.csv")
	write.csv(varR2, "Result/R2.csv")

# Constructing of allele matrix
#	if (Exp.Des == "diallel"){
#		allel.matrix <- Allele.matrix(mrk, selmrk, est.para, design = "diallel")
#	}
#	else if (Exp.Des == "NCII"){
#		allel.matrix <- Allele.matrix(mrk, selmrk, est.para, design = "NCII")
#	}
#	addmatrix <- allel.matrix$addMatrix
#	dommatrix <- allel.matrix$domMatrix

#	write.csv(addmatrix, "Result/Add_Matrix.csv")
#	write.csv(dommatrix, "Result/Dom_Matrix.csv")
# Prediction of the perforance of hybrid

}
