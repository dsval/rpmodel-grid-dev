#' Calculates the temperature dependence of the quantum yield efficiency
#'
#' Calculates the temperature dependence of the quantum yield efficiency
#' following the temperature dependence of the maximum quantum yield of photosystem II
#' in light-adapted tobacco leaves, determined by Bernacchi et al. (2003)
#'
#' @param tc Temperature, relevant for photosynthesis (degrees Celsius)
#' @param c4 Boolean specifying whether fitted temperature response for C4 plants
#' is used. Defaults to \code{FALSE} (C3 photoynthesis temperature resposne following
#' Bernacchi et al., 2003 is used).
#'
#' @details The temperature factor for C3 photosynthesis (argument \code{c4 = FALSE}) is calculated
#' based on Bernacchi et al. (2003) as
#' 			\deqn{
#' 				\phi(T) = 0.352 + 0.022 T - 0.00034 T^2
#'       }
#'
#' The temperature factor for C4 (argument \code{c4 = TRUE}) photosynthesis is calculated based on
#' unpublished work as
#' 			\deqn{
#' 				\phi(T) = -0.008 + 0.00375 T - 0.58e-4 T^2
#'       }
#'
#' The factor \eqn{\phi(T)} is to be multiplied with leaf absorptance and the fraction
#' of absorbed light that reaches photosystem II. In the P-model these additional factors
#' are lumped into a single apparent quantum yield efficiency parameter (argument \code{kphio}
#' to function \link{rpmodel}).
#'
#' @return A numeric value for \eqn{\phi(T)}
#'
#' @examples
#' ## Relative change in the quantum yield efficiency
#' ## between 5 and 25 degrees celsius (percent change):
#' print(paste((calc_ftemp_kphio(25.0)/calc_ftemp_kphio(5.0)-1)*100 ))
#'
#' @references  Bernacchi, C. J., Pimentel, C., and Long, S. P.:  In vivo temperature
#' 				response func-tions  of  parameters required  to  model  RuBP-limited
#' 				photosynthesis,  Plant  Cell Environ., 26, 1419–1430, 2003
#'
#' @export
#'
calc_ftemp_kphio <- function( tc, c4 = FALSE ){

	if (c4){
    ftemp = -0.008 + 0.00375 * tc - 0.58e-4 * tc^2   # Based on calibrated values by Shirley
	} else {
	  ftemp <- 0.352 + 0.022 * tc - 3.4e-4 * tc^2
	}

  ## avoid negative values
  ftemp <- max(ftemp, 0.0)

  return(ftemp)
}
