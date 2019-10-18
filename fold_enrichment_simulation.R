enrichment_calculation=function(coannotation){
  cumsum(coannotation)/ ( mean(coannotation)*(1:length(coannotation)) )
}



n=1000

x=c(rep(1,900),rep(0,n-900))

x1=c(rep(1,500),rep(0,n-500))

x2=c(rep(1,100),rep(0,n-100))

x3=c(rep(1,1),rep(0,n-1))

enrichment_x=enrichment_calculation(x)
enrichment_x1=enrichment_calculation(x1)
enrichment_x2=enrichment_calculation(x2)
enrichment_x3=enrichment_calculation(x3)

plot(1:length(x),enrichment_x)
max(enrichment_x)

plot(1:length(x1),enrichment_x1)
max(enrichment_x1)

plot(1:length(x2),enrichment_x2)
max(enrichment_x2)

plot(1:length(x3),enrichment_x3)
max(enrichment_x3)

#subtract random instead of dividing
enrichment_calculation_subtract=function(coannotation){
  cumsum(coannotation) - ( mean(coannotation)*(1:length(coannotation)) )
}


enrichment_subtract_x=enrichment_calculation_subtract(x)
enrichment_subtract_x1=enrichment_calculation_subtract(x1)
enrichment_subtract_x2=enrichment_calculation_subtract(x2)
enrichment_subtract_x3=enrichment_calculation_subtract(x3)

plot(1:length(x),enrichment_subtract_x)
max(enrichment_subtract_x)

plot(1:length(x1),enrichment_subtract_x1)
max(enrichment_subtract_x1)

plot(1:length(x2),enrichment_subtract_x2)
max(enrichment_subtract_x2)

plot(1:length(x3),enrichment_subtract_x3)
max(enrichment_subtract_x3)



