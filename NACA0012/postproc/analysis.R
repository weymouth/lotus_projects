library(dplyr);library(ggplot2)

ppdf = function(plot,name){
     pdf(name,8,4)
     print(plot)
     dev.off()
}
mod_fun <- function(E){Mod(fft(E))/length(E)*2}
k_fun <- function(z){
	l = length(z)
	dz = z[2]-z[1]
	2*pi/(dz*l)*(1:l-1)}

data = read.table("fort.10",col.names = c("time","z","E"))
data$E = data$E*75

sdat=data %>% group_by(time) %>%
		summarize(mean=mean(E),min=min(E),max=max(E))

q = ggplot(sdat,aes(time,mean,ymax=max,ymin=min))+
		geom_point()+geom_errorbar()+ylab(expression(E[3][D]))
ppdf(q,"04E_t.pdf")

freq=data %>% filter(time>mean(time)) %>% 
		group_by(time) %>%
		transmute(k=k_fun(z),mod=mod_fun(E)) %>%
		filter(k<max(k)/2,k>0,k<64)

q = ggplot(freq,aes(k,mod,group=k))+scale_y_log10()+
		geom_boxplot()+ylab(expression(tilde(E)[3][D]))
ppdf(q,"05E_k.pdf")
