#' NBER-CES Manufacturing Dairy Data
#'
#' Dataset of the National Bureau of Economic Research (NBER) and U.S. Census Bureau's Center for Economic Studies (CES). It contains information on the annual industry-level from 1958-2016 of the US.
#'
#' @format 
#' manuf is a data frame with 6188 rows and 7 columns:
#' \describe{
#'   \item{naics}{NAICS 2012 6-digit industry code}
#'   \item{year}{Year from 2000 to 2016}
#'   \item{Y}{Total value of shipments in millions of 2012 dollars}
#'   \item{K}{Real capital stock in millions of 2012 dollars}
#'   \item{L}{Production worker hours in millions}
#'   \item{M}{Total cost of materials  in millions of 2012 dollars}
#'   \item{I}{New capital spending in millions of 2012 dollars}
#' }
#' 
#' @details 
#' This is a subset of the data which contains data from 2000-2016 on output, employment, materials, investment and capital stocks.
#' 
#' @source <https://www.nber.org/research/data/nber-ces-manufacturing-industry-database>
#' 
#' @references
#' \insertRef{bartlesman1996nber}{dsfa}
#' 
"manuf"
