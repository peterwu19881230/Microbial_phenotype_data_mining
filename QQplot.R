##Understading QQplot

set.seed(0)
x=rnorm(10); print(x)
qqnorm(x)
qqnorm(2*x)
qqnorm(x+3)
qqline(x+3)
##The values of the normal distribution that is used here
c(qnorm(1/11),qnorm(2/11),qnorm(3/11),qnorm(4/11),qnorm(5/11),
  qnorm(6/11),qnorm(7/11),qnorm(8/11),qnorm(9/11),qnorm(10/11)
  )


n=rnorm(3)
qqnorm(n)
c(qnorm(1/4),qnorm(2/4),qnorm(3/4))
