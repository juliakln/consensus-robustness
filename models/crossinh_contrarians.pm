ctmc

// when volume changes, one should change all the parameters so to scale, including the percentage of zealots
const int N = 100; // 100
const int Nhalf = 50; // 50

const int m=50; //50
const int d=10; //10
const int dHalf=0;  // =5 for holding time

const int t=35;
const int h=40;


const int cHalf; // number of contrarians divided by two 
const int cAll = cHalf*2; //number of all contrarians

const double qx;
const double qy;

module cross_inhibition
	x : [0..N] init Nhalf-cHalf+dHalf;//x : [0..N] init Nhalf-cHalf;
	y : [0..N] init Nhalf-cHalf-dHalf;//y : [0..N] init Nhalf-cHalf;
	u : [0..N] init 0;//u : [0..N] init 0;
	Cx: [0..cAll] init cHalf;
	Cy: [0..cAll] init cHalf;

	[cix] 	   (x>0) & (y>0) & (u<N) -> x*y : (x'=x) & (y'=y-1) & (u'=u+1); // x+y -> x+u
	[ciy] 	   (x>0) & (y>0) & (u<N) -> x*y : (x'=x-1) & (y'=y) & (u'=u+1); // x+y -> y+u
	[rx] 	   (x>0) & (x<N) & (u>0) -> x*u : (x'=x+1) & (u'=u-1);  	// u+x -> 2x 
	[ry]       (y>0) & (y<N) & (u>0) -> y*u : (y'=y+1) & (u'=u-1);		// u+y -> 2y

 	[conxa]    (x>0) & (Cy>0) & (u<N) -> x*Cy : (x'=x-1) & (u'=u+1);		// x+Cy -> u+Cy
        [conxb]    (u>0) & (Cy>0) & (y<N) -> u*Cy : (u'=u-1) & (y'=y+1);		// u+Cy -> y+Cy
	[conxc]    (x>0) & (Cx>0) & (Cy<cAll) -> x*Cx : (Cx'=Cx-1) & (Cy'=Cy+1);	// x+Cx -> x+Cy
        [conya]    (y>0) & (Cx>0) & (u<N) -> y*Cx : (y'=y-1) & (u'=u+1);		// y+Cx -> u+Cx  (Cx<cAll) &
        [conyb]    (u>0) & (Cx>0) & (x<N) -> u*Cx : (u'=u-1) & (x'=x+1);		// u+Cx -> x+Cx
        [conyc]    (y>0) & (Cy>0) & (Cx<cAll) -> y*Cy : (Cy'=Cy-1) & (Cx'=Cx+1);	// y+Cy -> y+Cx
        [conxx]    (Cx>2) & (Cy<(cAll-1)) -> (Cx*(Cx-1)/2) : (Cx'=Cx-2) & (Cy'=Cy+2);	// Cx+Cx->Cy+Cy
        [conyy]    (Cy>2) & (Cx<(cAll-1)) -> (Cy*(Cy-1)/2) : (Cy'=Cy-2) & (Cx'=Cx+2);	// Cy+Cy->Cx+Cx


endmodule

// module representing the base rates of reactions
module base_rates
	
	[cix] true -> qx : true;
	[ciy] true -> qy : true;
	[rx] true -> qx : true;
	[ry] true -> qy : true;

        [conxa] true -> qy : true;
        [conxb] true -> qy : true;	
        [conxc] true -> qx : true;
        [conya] true -> qx : true;
        [conyb] true -> qx : true;
        [conyc] true -> qy : true;
        [conxx] true -> qx : true;
        [conyy] true -> qy : true;
endmodule


rewards "reaching_time"

//        (x+Cx>=m)&(x+Cx>=Cy+y+d): 1;//count time only if in xwinning states
//	!((x+Zx>=m)&(x>=y+d)): 1;
         !(((x+Cx>=m)&(x+Cx>=y+Cy+d)) | ((y+Cy>=m)&(y+Cy>=x+Cx+d))) : 1; // expected time to reach majority -> reward if no majority -- now without ! it's a reward for having majority
endrewards

rewards "holding_time"
	(((x+Cx>=m)&(x+Cx>=y+Cy+d)) | ((y+Cy>=m)&(y+Cy>=x+Cx+d))) : 1; 
endrewards


label "x_wins"=((x+Cx>=m)&(x+Cx>=y+Cy+d));
label "y_wins"=((y+Cy>=m)&(y+Cy>=x+Cx+d));
