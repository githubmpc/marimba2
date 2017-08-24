load_all()
mprob <- gMendelian()
mygsub <- function(x){
  x <- gssub("F_", "", x)
  x <- gsub("AA", "2", x)
  x <- gsub("AB", "1", x)
  x <- gsub("BB", "0", x)
  x
}
g1 <- mprob[, , 1] %>% as.tibble %>%
  gather(key="father_z", value="p") %>%
  mutate(father_z=mygsub(father_z),
         mother_z="0")
g2 <- mprob[, , 2] %>% as.tibble %>%
  gather(key="father_z", value="p") %>%
  mutate(father_z=mygsub(father_z),
         mother_z="1")
g3 <- mprob[, , 3] %>% as.tibble %>%
  gather(key="father_z", value="p") %>%
  mutate(father_z=mygsub(father_z),
         mother_z="2")
mprob <- bind_rows(g1, g2, g3) %>%
  set_colnames(c("f", "p", "m")) %>%
  select(c("f", "m", "p")) %>%
  mutate(parents=paste0(f, m),
         child_prob=rep(paste0("p(", 0:2, "|f,m)"), 9))
mprob2 <- mprob %>% select(c("parents", "child_prob", "p")) %>%
  spread(child_prob, p)
saveRDS(mprob2, file="../inst/extdata/mendelian_probs.rds")
