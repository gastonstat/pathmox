pathmox
============================

**pathmox** is the R package dedicated to the PATHMOX approach for obtaining segmentation trees in Partial Least Squares Path Modeling (PLS-PM). 


## Installation

Stable version on [CRAN](http://cran.r-project.org/web/packages/pathmox/index.html)
```ruby
install.packages("pathmox")
```

Development version on [github](https://github.com/gastonstat/pathmox)
```ruby
# install.packages("devtools") 
library(devtools)
install_github('pathmox', username='gastonstat')
```

## Example Usage with a Customer Satisfaction Model 
```ruby
library(pathmox)

# load data csimobile
data(csimobile)

# select manifest variables
data_mobile = csimobile[,8:33]
  
# path matrix (inner model)
IMAG = c(0, 0, 0, 0, 0, 0, 0)
EXPE = c(1, 0, 0, 0, 0, 0, 0)
QUAL = c(0, 1, 0, 0, 0, 0, 0)
VAL = c(0, 1, 1, 0, 0, 0, 0)
SAT = c(1, 1, 1, 1, 0, 0, 0)
COM = c(0, 0, 0, 0, 1, 0, 0)
LOY = c(1, 0, 0, 0, 1, 1, 0)
mob_path = rbind(IMAG, EXPE, QUAL, VAL, SAT, COM, LOY)

# list of blocks (outer model)
mob_blocks = list(1:5, 6:9, 10:15, 16:18, 19:21, 22:24, 25:26)

# vector of modes (reflective indicators)
mob_modes = rep("A", 7)

# apply plspm
mob_pls = plspm(data_mobile, mob_path, mob_blocks, modes = mob_modes, scheme = "factor", scaled = FALSE)

# plot inner model results
plot(mob_pls)

# re-ordering those segmentation variables with ordinal scale (Age and Education)
csimobile$Education = factor(csimobile$Education, 
    levels=c("basic","highschool","university"), ordered=TRUE)

# select the segmentation variables
seg_vars = csimobile[,1:7]

# Pathmox Analysis
mob_pathmox = pathmox(mob_pls, seg_vars, signif=.10, size=.10, deep=2)
```

More info at [www.gastonsanchez.com](http://www.gastonsanchez.com)

Links
-----
[pathmox package github](http://github.com/gastonstat/pathmox)


Author Contact
--------------
Gaston Sanchez (gaston.stat at gmail.com)
