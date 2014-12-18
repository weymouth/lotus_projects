require(ggplot2)

t = seq(0,1,0.01)
xp = 0.16
etap = 1-xp

U = function(t){1-t}
dU = function(xi){-xi}

alpha = function(t){pi/2.*(t-sin(2.*pi*t)/(2.*pi))}
dalpha = function(t,xi){pi/2.*(1-cos(2.*pi*t))*xi}
ddalpha = function(t,xi){pi^2*sin(2.*pi*t)*xi^2}

Uperp = function(t,xi){U(t)*sin(alpha(t))+etap*dalpha(t,xi)}
Upar = function(t){U(t)*cos(alpha(t))}
Utip = function(t,xi){sqrt(Uperp(t,xi)^2+Upar(t)^2)}
Uf = function(t,xi){2*U(t)*sin(alpha(t))+.55*dalpha(t,xi)}

dGamma = function(t,xi){Utip(t,xi)*Uperp(t,xi)}
Fadded = function(t,xi){pi/4*(dalpha(t,xi)*cos(alpha(t))*U(t)+sin(alpha(t))*dU(xi)+ddalpha(t,xi)*(.5-xp))}

Clg = function(t,xi){2.*dGamma(t,xi)*cos(alpha(t))}
Clp = function(t,xi){2.*Fadded(t,xi)*cos(alpha(t))}

Cl = function(t,xi){Clp(t,xi)+Clg(t,xi)}

Cdg = function(t,xi){2.*dGamma(t,xi)*sin(alpha(t))}
Cdp = function(t,xi){2.*Fadded(t,xi)*sin(alpha(t))}

Cd = function(t,xi){Cdp(t,xi)+Cdg(t,xi)}

