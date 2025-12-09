ctmc

// when volume changes, one should change all the parameters so to scale, including the percentage of zealots
const int N = 100;
const int Nhalf = 50;

const int m=50;
const int d=10;
const int t=35;
const int h=40;

const int Zx;
const int Zy;

const double qx; //1.05/N
const double qy; //0.95/N

module cross_inhibition
	x : [0..N] init Nhalf-Zx;
	y : [0..N] init Nhalf-Zy;
	u : [0..N] init 0;
	
	[cix] 	   (x>0) & (y>0) & (u<N) -> x*y : (x'=x) & (y'=y-1) & (u'=u+1); // x+y -> x+u
	[ciy] 	   (x>0) & (y>0) & (u<N) -> x*y : (x'=x-1) & (y'=y) & (u'=u+1); // x+y -> y+u
	[rx] 	   (x>0) & (x<N) & (u>0) -> x*u : (x'=x+1) & (u'=u-1);  	// u+x -> 2x 
	[ry]       (y>0) & (y<N) & (u>0) -> y*u : (y'=y+1) & (u'=u-1);		// u+y -> 2y

	[zeaxa]    (y>0) & (u<N) & (Zx>0) -> y*Zx : (y'=y-1) & (u'=u+1);		// y+Zx -> u+Zx
	[zeaxb]    (u>0) & (x<N) & (Zx>0) -> u*Zx : (x'=x+1) & (u'=u-1);		// u+Zx -> x+Zx
	[zeaya]    (x>0) & (u<N) & (Zy>0) -> x*Zy : (x'=x-1) & (u'=u+1);		// x+Zy -> u+Zy
	[zeayb]    (u>0) & (y<N) & (Zy>0) -> u*Zy : (y'=y+1) & (u'=u-1);		// u+Zy -> y+Zy

endmodule

// module representing the base rates of reactions
module base_rates
	
	[cix] true -> qx : true;
	[ciy] true -> qy : true;
	[rx] true -> qx : true;
	[ry] true -> qy : true;

	[zeaxa] true -> qx : true;	
	[zeaxb] true -> qx : true;
	[zeaya] true -> qy : true;
	[zeayb] true -> qy : true;

endmodule

rewards "reaching_time"

//        (x+Cx>=m)&(x+Cx>=Cy+y+d): 1;//count time only if in xwinning states
//	!((x+Zx>=m)&(x>=y+d)): 1;
         !(((x+Zx>=m)&(x>=y+d)) | ((y+Zy>=m)&(y>=x+d))) : 1; // expected time to reach majority -> reward if no majority -- now without ! it's a reward for having majority
endrewards

rewards "holding_time"
	(((x+Zx>=m)&(x>=y+d)) | ((y+Zy>=m)&(y>=x+d))) : 1;
endrewards



label "x_wins"=((x+Zx>=m)&(x>=y+d));
label "y_wins"=((y+Zy>=m)&(y>=x+d));
