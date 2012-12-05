
famax=@(x) sqrt(2)*(sin(x)-1)+cos(x)
alphamax=180*(1-1/pi*fzero(famax,30*pi/180,optimset('Display','iter')))