NOIA_A <- function(Code)
{
  fg=table(Code)/length(Code)
  Xa=ifelse(Code == 0, fg[2]+2*fg[3], ifelse(Code == 1, -(1-fg[2]-2*fg[3]), -(2-fg[2]-2*fg[3])))
  return(Xa)
}

NOIA_D <- function(Code)
{
  fg=table(Code)/length(Code)
  Xd=ifelse(Code == 0, -2*fg[2]*fg[3], ifelse(Code == 1, 4*fg[3]*fg[1], -2*fg[2]*fg[1]))
  return(Xd)
}
