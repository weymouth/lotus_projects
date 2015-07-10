library(dplyr);library(ggplot2)

ppdf = function(plot,name){
     pdf(name,8,4)
     print(plot)
     dev.off()
}

mod_fun <- function(E){Mod(fft(E))/length(E)*2}
k_fun <- function(z){2*pi/z[length(z)]*seq(0,length(z)-1)}

read_dat <-function(){
    tab = read.table("fort.11",col.names = c("time","x","y","z","u","v","w"))

    tab %>% mutate(ry = round(y,2), rt = round(time)) %>%
		group_by(y,time) %>%
		mutate(e=0.5*((u-mean(u))^2+(v-mean(v))^2+(w-mean(w))^2))
}

data = read_dat()

#q = ggplot(data,aes(u,y,group=y))+geom_jitter(alpha=0.3)+facet_wrap( ~ rt,nrow=3)
#ppdf(q,'05_u_y_t.pdf')

#q = ggplot(data,aes(e,y,group=y))+geom_jitter(alpha=0.3)+facet_wrap( ~ rt,nrow=3)
#ppdf(q,'06_e_y_t.pdf')

q = qplot(round(time),e,group=time,data=subset(data,y>0),geom='boxplot',log='y')+
	labs(x='time',y='E')
ppdf(q,"04e_t.pdf")

freq=data %>% subset(time>mean(time),y>0) %>%
		group_by(y,time) %>%
		transmute(k=k_fun(z),E=mod_fun(e)) %>%
		filter(k<max(k)/2,k>0)

q = qplot(k,E,data=freq,log='xy',geom='smooth')+ylab('')+
  geom_smooth(aes(y=k*k*E*2/5.3e3),color='red')
ppdf(q,"05E_k.pdf")

library(animation)

FUN <- function(t){
	ggplot(subset(data,rt==t),aes(z,u,color=y,group=y))+
	geom_line()+ylim(range(data$u))+theme_minimal()
}

#saveGIF({for(i in unique(data$rt)) print(FUN(i))},ani.width=960)
