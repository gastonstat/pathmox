\name{csimobile}
\alias{csimobile}
\docType{data}
\title{CSI Mobile Phone dataset}
\description{
  This table contains data from a Customer Satisfaction Index CSI on Mobile Phone
}
\usage{
data(csimobile)
}
\format{
  A data frame with 87 observations on 33 variables. The first seven variables are segmentation variables. The rest of the variables refer to seven latent concepts: 1) \code{IMAG}=Image, 2) \code{EXPE}=Expectations, 3) \code{QUAL}=Quality, 4) \code{VAL}=Value, 5) \code{SAT}=Satisfaction, 6) \code{COM}=Complaints, and 7) \code{LOY}=Loyalty. 

  \code{IMAG}: Includes variables such as trustworthiness, dynamic, solidness, innovation, and caring about customer's needs. \cr
  \code{EXPE}: Includes variables such as products and services provided and expectations for the overall quality. \cr
  \code{QUAL}: Includes variables such as reliable products and services, range of products and services, and overall perceived quality. \cr
  \code{VAL}:  Includes variables such as service and products, quality relative to price, and price relative to quality. \cr
  \code{SAT}:  Includes variables such as overall rating of satisfaction, fulfillment of expectations, satisfaction relative to other phone providers.\cr
  \code{COM}:  Includes one variable defining how well or poorly custmer's complaints were handled. \cr
  \code{LOY}:  Includes variables such as propensity to choose the same phone provider again, intention to recommend the phone provider to friends. \cr
  
  \describe{
    \item{\code{Age}}{a factor with levels \code{<=25} \code{>25}}
    \item{\code{Gender}}{a factor with levels \code{female} \code{male}}
    \item{\code{Occupation}}{a factor with levels \code{employee} \code{other} \code{student}}
    \item{\code{Education}}{a factor with levels \code{basic} \code{highschool} \code{university}}
    \item{\code{Operator}}{a factor with levels \code{A} \code{B} \code{C}}
    \item{\code{Contract}}{a factor with levels \code{contract} \code{prepay}}
    \item{\code{SwitchOp}}{a factor with levels \code{No } \code{Yes }}
    \item{\code{image1}}{First MV of the block Image}
    \item{\code{image2}}{Second MV of the block Image}
    \item{\code{image3}}{Thrid MV of the block Image}
    \item{\code{image4}}{Fourth MV of the block Image}
    \item{\code{image5}}{Fifth MV of the block Image}
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
    \item{\code{value1}}{First MV of the block Value}
    \item{\code{value2}}{Second MV of the block Value}
    \item{\code{value3}}{Third MV of the block Value}
    \item{\code{satis1}}{First MV of the block Satisfaction}
    \item{\code{satis2}}{Second MV of the block Satisfaction}
    \item{\code{satis3}}{Third MV of the block Satisfaction}
    \item{\code{comp1}}{First MV of the block Complaints}
    \item{\code{comp2}}{Second MV of the block Complaints}
    \item{\code{comp3}}{Third MV of the block Complaints}
    \item{\code{loyal1}}{First MV of the block Loyalty}
    \item{\code{loyal2}}{Second MV of the block Loyalty}
  }
}
\source{
  Laboratory of Information Analysis and Modeling (LIAM). Facultat d'Informatica de Barcelona, Universitat Politecnica de Catalunya.
}
\examples{
data(csimobile)
}
\keyword{datasets}
