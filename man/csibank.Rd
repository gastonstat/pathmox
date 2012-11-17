\name{csibank}
\alias{csibank}
\docType{data}
\title{CSI Bank dataset}
\description{
  This data set contains the variables from a customer satisfaction study of a Spanish credit institution.
}
\usage{
data(csibank)
}
\format{
  A data frame with 1707 observations on the following 32 variables. The first five variables are segmentation variables. The rest of the variables refer to seven latent concepts: 1) \code{IMAG}=Image, 2) \code{EXPE}=Expectations, 3) \code{QUAL}=Quality, 4) \code{VAL}=Value, 5) \code{SAT}=Satisfaction, and 6) \code{LOY}=Loyalty. 

    \code{IMAG}: Includes variables such as reputation, trustworthiness, seriousness, solidness, and caring about customer's needs. \cr
    \code{EXPE}: Includes variables such as products and services provided, customer service, providing solutions, and expectations for the overall quality. \cr
    \code{QUAL}: Includes variables such as reliable products and services, range of products and services, personal advice, and overall perceived quality. \cr
    \code{VAL}: Includes variables such as beneficial services and products, valuable investments, quality relative to price, and price relative to quality. \cr
    \code{SAT}: Includes variables such as overall rating of satisfaction, fulfillment of expectations, satisfaction relative to other banks, and performance relative to customer's ideal bank. \cr
    \code{LOY}: Includes variables such as propensity to choose the same bank again, propensity to switch to other bank, intention to recommend the bank to friends, and sense of loyalty. \cr
  
  \describe{
    \item{\code{Gender}}{a factor with levels \code{Female} \code{Male}}
    \item{\code{Age}}{a factor with levels \code{<=25} \code{>=66} \code{26-35} \code{36-45} \code{46-55} \code{56-65}}
    \item{\code{Education}}{a factor with levels \code{Elementary} \code{Graduated} \code{Highschool} \code{Undergrad} \code{Unfinished}}
    \item{\code{Occupation}}{a factor with levels \code{Manager} \code{MediumEmplo} \code{Notemploy} \code{OwnFreelan} \code{Retired}}
    \item{\code{Region}}{a factor with levels \code{Center} \code{East} \code{North}}
    \item{\code{imag1}}{First MV of the block Image}
    \item{\code{imag2}}{Second MV of the block Image}
    \item{\code{imag3}}{Third MV of the block Image}
    \item{\code{imag4}}{Fourth MV of the block Image}
    \item{\code{imag5}}{Fifth MV of the block Image}
    \item{\code{imag6}}{Sixth MV of the block Image}
    \item{\code{expe1}}{First MV of the block Expectations}
    \item{\code{expe2}}{Second MV of the block Expectations}
    \item{\code{expe3}}{Third MV of the block Expectations}
    \item{\code{expe4}}{Fourth MV of the block Expectations}
    \item{\code{qual1}}{First MV of the block Quality}
    \item{\code{qual2}}{Second MV of the block Quality}
    \item{\code{qual3}}{Third MV of the block Quality}
    \item{\code{qual4}}{Fourth MV of the block Quality}
    \item{\code{qual5}}{Fifth MV of the block Quality}
    \item{\code{qual6}}{Sixth MV of the block Quality}
    \item{\code{qual7}}{Seventh MV of the block Quality}
    \item{\code{val1}}{First MV of the block Value}
    \item{\code{val2}}{Second MV of the block Value}
    \item{\code{val3}}{Third MV of the block Value}
    \item{\code{val4}}{Fourth MV of the block Value}
    \item{\code{sat1}}{First MV of the block Satisfaction}
    \item{\code{sat2}}{Second MV of the block Satisfaction}
    \item{\code{sat3}}{Third MV of the block Satisfaction}
    \item{\code{loy1}}{First MV of the block Loyalty}
    \item{\code{loy2}}{Second MV of the block Loyalty}
    \item{\code{loy3}}{Third MV of the block Loyalty}
  }
}
\source{
  Laboratory of Information Analysis and Modeling (LIAM). Facultat d'Informatica de Barcelona, Universitat Politecnica de Catalunya.
}
\examples{
data(csibank)
}
\keyword{datasets}
