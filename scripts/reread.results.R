#re-read in old results for additional analyses

theta.tbl <- read.table("theta.10.1000.c1.txt", sep="\t")
sigma.tbl <- read.table("sigma.10.1000.c1.txt", sep="\t")
mixture.tbl <- read.table("mixture.parents.10.1000.c1.txt", sep="\t")
mixture.tbl.offspring <- read.table("mixture.offspring.10.1000.c1.txt", sep="\t")
parent.misclass.tbl <- read.table("misclas.par.10.1000.c1.txt", sep="\t")
child.misclass.tbl <- read.table("misclas.child.10.1000.c1.txt", sep="\t")
misclass.perc <- read.table("misclas.perc.10.1000.c1.txt", sep="\t")
prop.mis.tbl <- read.table("misclass.proportion.ci.c1.txt", sep="\t")
prop.mis.tbl.q <- read.table("misclass.proportion.ci.quantile.c1.txt", sep="\t")
mistake.typed.tbl <- read.table("mistake.typed.tbl.c1.txt", sep="\t")
ess.tbl <- read.table("ess.tbl.c1.txt", sep="\t")
mistake.values.tbl <- read.table("mistake.values.tbl.c1.txt", sep="\t")
