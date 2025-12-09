"""
Define voter model, cross-inhibition model, combined model, SC and CC (Cardelli paper) as RML models for Plasmalab.
Rates qx and qy are scaled by population size N.

Save in '../models/consensus_model.rml"

"""

import os


# change to location of this script
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

# define file for model 
model = "../models/consensus_model.rml"


"""
VOTER MODEL
    stubborn: z (zealots), c (contrarians)
    N: total population size
    init: initial number of stubborns for each opinion
    ratex: rate for opinion x
    ratey: rate for opinion y

    Output: RML file model
"""
def write_voter_model(stubborn, N, init, ratex = 1.0, ratey = 1.0):
    # compute remaining number of pure agents
    X = int(1/2 * (N - init))
    Z = int(1/2 * init)
    f = open(model, "w")
    if stubborn == 'z':
        f.write("""ctmc

            const int Zx = """ + str(Z) + """;
            const int Zy = """ + str(Z) + """;
            const int N = """ + str(N) + """;

            module cross_inhibition
                
                x : [0..N] init """ + str(X) + """;
                y : [0..N] init """ + str(X) + """;
                
                [switch1] 	 (x>0) & (y>0) & (x<N)  -> x*y : (x'=x+1) & (y'=y-1);       // x+y -> x+x 
                [switch2] 	 (x>0) & (y>0) & (y<N)  -> x*y : (x'=x-1) & (y'=y+1);       // x+y -> y+y
                [zealotx]    (y>0) & (x<N) & (Zx>0) -> y*Zx : (y'=y-1) & (x'=x+1);		// y+Zx -> x+Zx
                [zealoty]    (x>0) & (y<N) & (Zy>0)	-> x*Zy : (x'=x-1) & (y'=y+1);		// x+Zy -> y+Zy
                [dead1]      (x=0) -> (x'=x);
                [dead2]      (y=0) -> (y'=y);

            endmodule

            // base rates
            const double qx = """ + str(ratex/N) + """; 
            const double qy = """ + str(ratey/N) + """; 

            // module representing the base rates of reactions
            module base_rates
                
                [switch1] true -> qx : true;
                [switch2] true -> qy : true;
                [zealotx] true -> qx : true;	
                [zealoty] true -> qy : true;

            endmodule""")
    elif stubborn == 'c':
        f.write("""ctmc

            const int N = """ + str(N) + """;
            const int cN = """ + str(init) + """;

            module cross_inhibition
                
                x : [0..N] init """ + str(X) + """;
                y : [0..N] init """ + str(X) + """;
                Cx : [0..N] init """ + str(Z) + """;
                Cy : [0..N] init """ + str(Z) + """;
                
                [cix] 	   (x>0) & (y>0) & (x<N)  -> x*y : (x'=x+1) & (y'=y-1);  // x+y -> x+x
                [ciy] 	   (x>0) & (y>0) & (y<N)  -> x*y : (x'=x-1) & (y'=y+1);  // x+y -> y+y

                [conxa]    (x>0) & (Cy>0) & (y<N) -> x*Cy : (x'=x-1) & (y'=y+1);		    // x+Cy -> y+Cy
                [conxc]    (x>0) & (Cx>0) & (Cy<N) -> x*Cx : (Cx'=Cx-1) & (Cy'=Cy+1);	    // x+Cx -> x+Cy
                [conya]    (y>0) & (Cx>0) & (x<N) -> y*Cx : (y'=y-1) & (x'=x+1);		    // y+Cx -> x+Cx
                [conyc]    (y>0) & (Cy>0) & (Cx<N) -> y*Cy : (Cy'=Cy-1) & (Cx'=Cx+1);	    // y+Cy -> y+Cx
                [conxx]    (Cx>1) & (Cy<(N-1)) -> (Cx*(Cx-1)/2) : (Cx'=Cx-2) & (Cy'=Cy+2); // Cx+Cx->Cy+Cy
                [conyy]    (Cy>1) & (Cx<(N-1)) -> (Cy*(Cy-1)/2) : (Cy'=Cy-2) & (Cx'=Cx+2);	// Cy+Cy->Cx+Cx

            endmodule

            // base rates
            const double qx = """ + str(ratex/N) + """; 
            const double qy = """ + str(ratey/N) + """; 

            // module representing the base rates of reactions
            module base_rates
                
                [cix] true -> qx : true;
                [ciy] true -> qy : true;
                [conxa] true -> qy : true;
                [conxc] true -> qx : true;
                [conya] true -> qx : true;
                [conyc] true -> qy : true;
                [conxx] true -> qx : true;
                [conyy] true -> qy : true;

            endmodule""")

    else:
        raise Exception('Type of stubborn individual not supported.')
    f.close()


"""
CROSS-INHIBITION MODEL
    stubborn: z (zealots), c (contrarians), both
    N: total population size
    init: initial number of stubborns for each opinion
    init2: initial number of contrarians for each opinion (only if stubborn = both)
    ratex: rate for opinion x
    ratey: rate for opinion y

    Output: RML file model
"""
def write_crossinh_model(stubborn, N, init, init2=0, ratex = 1.0, ratey = 1.0):
    # compute remaining number of pure agents
    X = int(1/2 * (N - init))
    Z = int(1/2 * init)
    f = open(model, "w")
    if stubborn == 'z':
        f.write("""ctmc

            const int Zx = """ + str(Z) + """;
            const int Zy = """ + str(Z) + """;
            const int N = """ + str(N) + """;

            module cross_inhibition
                
                x : [0..N] init """ + str(X) + """;
                y : [0..N] init """ + str(X) + """;
                u : [0..N] init 0;
                
                [cix] 	   (x>0) & (y>0) & (u<N) -> x*y : (x'=x) & (y'=y-1) & (u'=u+1); // x+y -> x+u
                [ciy] 	   (x>0) & (y>0) & (u<N) -> x*y : (x'=x-1) & (y'=y) & (u'=u+1); // x+y -> y+u
                [rx] 	   (x>0) & (x<N) & (u>0) -> x*u : (x'=x+1) & (u'=u-1);  	// u+x -> 2x 
                [ry]       (y>0) & (y<N) & (u>0) -> y*u : (y'=y+1) & (u'=u-1);		// u+y -> 2y
                [zeaxa]    (y>0) & (u<N)	 -> y*Zx : (y'=y-1) & (u'=u+1);		// y+Zx -> u+Zx
                [zeaxb]    (u>0) & (x<N)	 -> u*Zx : (x'=x+1) & (u'=u-1);		// u+Zx -> x+Zx
                [zeaya]    (x>0) & (u<N)	 -> x*Zy : (x'=x-1) & (u'=u+1);		// x+Zy -> u+Zy
                [zeayb]    (u>0) & (y<N)	 -> u*Zy : (y'=y+1) & (u'=u-1);		// u+Zy -> y+Zy

            endmodule

            // base rates
            const double qx = """ + str(ratex/N) + """; 
            const double qy = """ + str(ratey/N) + """; 

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

            endmodule""")
    elif stubborn == 'c':
        f.write("""ctmc

            const int N = """ + str(N) + """;
            const int cN = """ + str(init) + """;

            module cross_inhibition
                
                x : [0..N] init """ + str(X) + """;
                y : [0..N] init """ + str(X) + """;
                u : [0..N] init 0;
                Cx : [0..N] init """ + str(Z) + """;
                Cy : [0..N] init """ + str(Z) + """;
                
                [cix] 	   (x>0) & (y>0) & (u<N) -> x*y : (x'=x) & (y'=y-1) & (u'=u+1); // x+y -> x+u
                [ciy] 	   (x>0) & (y>0) & (u<N) -> x*y : (x'=x-1) & (y'=y) & (u'=u+1); // x+y -> y+u
                [rx] 	   (x>0) & (x<N) & (u>0) -> x*u : (x'=x+1) & (u'=u-1);  	// u+x -> 2x 
                [ry]       (y>0) & (y<N) & (u>0) -> y*u : (y'=y+1) & (u'=u-1);		// u+y -> 2y

                [conxa]    (x>0) & (Cy>0) & (u<N) -> x*Cy : (x'=x-1) & (u'=u+1);		// x+Cy -> u+Cy
                [conxb]    (u>0) & (Cy>0) & (y<N) -> u*Cy : (u'=u-1) & (y'=y+1);		// u+Cy -> y+Cy
                [conxc]    (x>0) & (Cx>0) & (Cy<cN) -> x*Cx : (Cx'=Cx-1) & (Cy'=Cy+1);		// x+Cx -> x+Cy
                [conya]    (y>0) & (Cx>0) & (u<N) -> y*Cx : (y'=y-1) & (u'=u+1);		// y+Cx -> u+Cx
                [conyb]    (u>0) & (Cx>0) & (x<N) -> u*Cx : (u'=u-1) & (x'=x+1);		// u+Cx -> x+Cx
                [conyc]    (y>0) & (Cy>0) & (Cx<cN) -> y*Cy : (Cy'=Cy-1) & (Cx'=Cx+1);		// y+Cy -> y+Cx
                [conxx]    (Cx>2) & (Cy<(cN-1)) -> (Cx*(Cx-1)/2) : (Cx'=Cx-2) & (Cy'=Cy+2);		// Cx+Cx->Cy+Cy
                [conyy]    (Cy>2) & (Cx<(cN-1)) -> (Cy*(Cy-1)/2) : (Cy'=Cy-2) & (Cx'=Cx+2);		// Cy+Cy->Cx+Cx

            endmodule

            // base rates
            const double qx = """ + str(ratex/N) + """; 
            const double qy = """ + str(ratey/N) + """; 

            // module representing the base rates of reactions
            module base_rates
                
                [cix] true -> qx : true;
                [ciy] true -> qy : true;
                [rx] true -> qx : true;
                [ry] true -> qy : true;
                [conxa] true -> qy : true;
                [conxb] true -> qy : true;	
                [conxc] true -> qy : true;
                [conya] true -> qx : true;
                [conyb] true -> qx : true;
                [conyc] true -> qx : true;
                [conxx] true -> qy : true;
                [conyy] true -> qx : true;

            endmodule""")
        
    elif stubborn == 'both':
        f.write("""ctmc

        const int Zx = """ + str(Z) + """;
        const int Zy = """ + str(Z) + """;
        const int N = """ + str(N) + """;
        const int cN = """ + str(init2) + """;

        module cross_inhibition
            
            x : [0..N] init """ + str(X) + """;
            y : [0..N] init """ + str(X) + """;
            u : [0..N] init 0;
            Cx : [0..N] init """ + str(int(1/2 * init2)) + """;
	        Cy : [0..N] init """ + str(int(1/2 * init2)) + """;
            
            [cix] 	   (x>0) & (y>0) & (u<N) -> x*y : (x'=x) & (y'=y-1) & (u'=u+1); // x+y -> x+u
            [ciy] 	   (x>0) & (y>0) & (u<N) -> x*y : (x'=x-1) & (y'=y) & (u'=u+1); // x+y -> y+u
            [rx] 	   (x>0) & (x<N) & (u>0) -> x*u : (x'=x+1) & (u'=u-1);  	// u+x -> 2x 
            [ry]       (y>0) & (y<N) & (u>0) -> y*u : (y'=y+1) & (u'=u-1);		// u+y -> 2y

            [zeaxa]    (y>0) & (u<N)	 -> y*Zx : (y'=y-1) & (u'=u+1);		// y+Zx -> u+Zx
            [zeaxb]    (u>0) & (x<N)	 -> u*Zx : (x'=x+1) & (u'=u-1);		// u+Zx -> x+Zx
            [zeaya]    (x>0) & (u<N)	 -> x*Zy : (x'=x-1) & (u'=u+1);		// x+Zy -> u+Zy
            [zeayb]    (u>0) & (y<N)	 -> u*Zy : (y'=y+1) & (u'=u-1);		// u+Zy -> y+Zy

            [conxa]    (x>0) & (Cy>0) & (u<N) -> x*Cy : (x'=x-1) & (u'=u+1);		// x+Cy -> u+Cy
            [conxb]    (u>0) & (Cy>0) & (y<N) -> u*Cy : (u'=u-1) & (y'=y+1);		// u+Cy -> y+Cy
            [conxc]    (x>0) & (Cx>0) & (Cy<cN) -> x*Cx : (Cx'=Cx-1) & (Cy'=Cy+1);		// x+Cx -> x+Cy
            [conya]    (y>0) & (Cx>0) & (u<N) -> y*Cx : (y'=y-1) & (u'=u+1);		// y+Cx -> u+Cx
            [conyb]    (u>0) & (Cx>0) & (x<N) -> u*Cx : (u'=u-1) & (x'=x+1);		// u+Cx -> x+Cx
            [conyc]    (y>0) & (Cy>0) & (Cx<cN) -> y*Cy : (Cy'=Cy-1) & (Cx'=Cx+1);		// y+Cy -> y+Cx
            [conxx]    (Cx>2) & (Cy<(cN-1)) -> Cx*(Cx-1)/2 : (Cx'=Cx-2) & (Cy'=Cy+2);		// Cx+Cx->Cy+Cy
            [conyy]    (Cy>2) & (Cx<(cN-1)) -> Cy*(Cy-1)/2 : (Cy'=Cy-2) & (Cx'=Cx+2);		// Cy+Cy->Cx+Cx

            [zeaconx]  (Cx>0) & (Cy<cN) -> Cx*Zx : (Cx'=Cx-1) & (Cy'=Cy+1);      // Cx+Zx -> Cy+Zx
            [zeacony]  (Cy>0) & (Cx<cN) -> Cy*Zy : (Cy'=Cy-1) & (Cx'=Cx+1);      // Cy+Zy -> Cx+Zy
            
        endmodule

        // base rates
            const double qx = """ + str(ratex/N) + """; 
            const double qy = """ + str(ratey/N) + """; 

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
            [conxa] true -> qy : true;
            [conxb] true -> qy : true;	
            [conxc] true -> qy : true;
            [conya] true -> qx : true;
            [conyb] true -> qx : true;
            [conyc] true -> qx : true;
            [conxx] true -> qy : true;
            *[conyy] true -> qx : true;
            [zeaconx] true -> qy : true;
            [zeacony] true -> qx : true;

        endmodule""")

    else:
        raise Exception('Type of stubborn individual not supported.')
    f.close()



"""
COMBINED MODEL
-> contains 5x direct switch reactions from voter model & 1x indirect switch reactions from cross-inhibition model
"""
def write_combined_model(stubborn, N, init, ratex = 1.0, ratey = 1.0):
    # compute remaining number of pure agents
    X = int(1/2 * (N - init))
    Z = int(1/2 * init)
    f = open(model, "w")
    if stubborn == 'z':
        f.write("""ctmc

            const int Zx = """ + str(Z) + """;
            const int Zy = """ + str(Z) + """;
            const int N = """ + str(N) + """;

            module cross_inhibition
                
                x : [0..N] init """ + str(X) + """;
                y : [0..N] init """ + str(X) + """;
                u : [0..N] init 0;
                
                [voterx]   (x>0) & (y>0) & (x<N)  -> x*y : (x'=x+1) & (y'=y-1);  // x+y -> x+x
                [votery]   (x>0) & (y>0) & (y<N)  -> x*y : (x'=x-1) & (y'=y+1);  // x+y -> y+y
                [voterx2]   (x>0) & (y>0) & (x<N)  -> x*y : (x'=x+1) & (y'=y-1);  // x+y -> x+x
                [votery2]   (x>0) & (y>0) & (y<N)  -> x*y : (x'=x-1) & (y'=y+1);  // x+y -> y+y
                [voterx3]   (x>0) & (y>0) & (x<N)  -> x*y : (x'=x+1) & (y'=y-1);  // x+y -> x+x
                [votery3]   (x>0) & (y>0) & (y<N)  -> x*y : (x'=x-1) & (y'=y+1);  // x+y -> y+y
                [voterx4]   (x>0) & (y>0) & (x<N)  -> x*y : (x'=x+1) & (y'=y-1);  // x+y -> x+x
                [votery4]   (x>0) & (y>0) & (y<N)  -> x*y : (x'=x-1) & (y'=y+1);  // x+y -> y+y
                [voterx5]   (x>0) & (y>0) & (x<N)  -> x*y : (x'=x+1) & (y'=y-1);  // x+y -> x+x
                [votery5]   (x>0) & (y>0) & (y<N)  -> x*y : (x'=x-1) & (y'=y+1);  // x+y -> y+y
                 
                [cix] 	   (x>0) & (y>0) & (u<N) -> x*y : (x'=x) & (y'=y-1) & (u'=u+1); // x+y -> x+u
                [ciy] 	   (x>0) & (y>0) & (u<N) -> x*y : (x'=x-1) & (y'=y) & (u'=u+1); // x+y -> y+u

                [rx] 	   (x>0) & (x<N) & (u>0) -> x*u : (x'=x+1) & (u'=u-1);  	// u+x -> 2x 
                [ry]       (y>0) & (y<N) & (u>0) -> y*u : (y'=y+1) & (u'=u-1);	// u+y -> 2y

                [deadx] (x=0) -> (x'=x);                
                [deady] (y=0) -> (y'=y);

            endmodule

            // base rates
            const double qx = """ + str(ratex/N) + """; 
            const double qy = """ + str(ratey/N) + """; 
            const double cx = """ + str(ratex/N) + """; 
            const double cy = """ + str(ratey/N) + """; 

            // module representing the base rates of reactions
            module base_rates
                
                [voterx] true -> qx : true;
                [voterx2] true -> qx : true;
                [voterx3] true -> qx : true;
                [voterx4] true -> qx : true;
                [voterx5] true -> qx : true;
                [votery] true -> qy : true;
                [votery2] true -> qy : true;
                [votery3] true -> qy : true;
                [votery4] true -> qy : true;
                [votery5] true -> qy : true;
                [cix] true -> cx : true;
                [ciy] true -> cy : true;
                [rx] true -> cx : true;
                [ry] true -> cy : true;

            endmodule""")
    else:
        raise Exception('Type of stubborn individual not supported.')
    f.close()



"""
SC model taken from Cardelli paper
    N: total population size
    ratex: rate for opinion x
    ratey: rate for opinion y
    ratet: rate for all other reactions

    Output: RML file model
"""
def write_sc_model(N, ratex = 1.0, ratey = 1.0, ratet = 1.0):
    # compute remaining number of pure agents
    X = int(1/3 * N)
    f = open(model, "w")
    f.write("""ctmc

        const int N = """ + str(N) + """;
        const int Nhalf = """ + str(X) + """;

        const int xMax = xInit+yInit+bInit;
        const int yMax = xInit+yInit+bInit;
        const int bMax = xInit+yInit+bInit;
        const int wMax = wInit+zInit+uInit;
        const int zMax = wInit+zInit+uInit;
        const int uMax = wInit+zInit+uInit;
        const int pMax = pInit+rInit+qInit;
        const int rMax = pInit+rInit+qInit;
        const int qMax = pInit+rInit+qInit;
        const int sMax = wInit+zInit+uInit;
        const int tMax = pInit+rInit+qInit;

        const int xInit = Nhalf;
        const int yInit = Nhalf;
        const int bInit = Nhalf;
        const int wInit = Nhalf;
        const int zInit = Nhalf;
        const int uInit = Nhalf;
        const int pInit = Nhalf;
        const int rInit = Nhalf;
        const int qInit = Nhalf;
        const int sInit = Nhalf;
        const int tInit = Nhalf;

        module SC

            x : [0..xMax] init xInit;
            y : [0..yMax] init yInit;
            b : [0..bMax] init bInit;
            w : [0..wMax] init wInit;
            z : [0..zMax] init zInit;
            u : [0..uMax] init uInit;
            p : [0..pMax] init pInit;
            r : [0..rMax] init rInit;
            q : [0..qMax] init qInit;
            s : [0..sMax] init sInit;
            t : [0..tMax] init tInit;

            // X + W -> W + B
            [sc1] x>0 & w>0 -> (x*w) : (x'=x-1) & (b'=min(b+1, bMax));
            // B + W -> W + Y
            [sc2] b>0 & w>0 -> (b*w) : (b'=b-1) & (y'=min(y+1, yMax));
            // Y + R -> R + B
            [sc3] y>0 & r>0 -> (y*r) : (y'=y-1) & (b'=min(b+1, bMax));
            // B + R -> R + X
            [sc4] b>0 & r>0 -> (b*r) : (b'=b-1) & (x'=min(x+1, xMax));

            // W + S -> S + U
            [sc5] w>0 & s>0 -> (w*s) : (w'=w-1) & (u'=min(u+1, uMax));
            // U + S -> S + Z
            [sc6] u>0 & s>0 -> (u*s) : (u'=u-1) & (z'=min(z+1, zMax));
            // Z + Y -> Y + U
            [sc7] z>0 & y>0 -> (z*y) : (z'=z-1) & (u'=min(u+1, uMax));
            // U + Y -> Y + W
            [sc8] u>0 & y>0 -> (u*y) : (u'=u-1) & (w'=min(w+1, wMax));

            // P + X -> X + Q
            [sc9] p>0 & x>0 -> (p*x) : (p'=p-1) & (q'=min(q+1, qMax));
            // Q + X -> X + R
            [sca1] q>0 & x>0 -> (q*x) : (q'=q-1) & (r'=min(r+1, rMax));
            // R + T -> T + Q
            [sca2] r>0 & t>0 -> (r*t) : (r'=r-1) & (q'=min(q+1, qMax));
            // Q + T -> T + P
            [sca3] q>0 & t>0 -> (q*t) : (q'=q-1) & (p'=min(p+1, pMax));

        endmodule

        // base rates
        const double rx = """ + str(ratex/N) + """;
        const double ry = """ + str(ratey/N) + """;
        const double rt = """ + str(ratet/N) + """;

        // module representing the base rates of reactions
        module base_rates
            
        [sc1] true -> ry: true;
        [sc2] true -> ry: true;
        [sc3] true -> rx: true;
        [sc4] true -> rx: true;
        [sc5] true -> rt: true;
        [sc6] true -> rt: true;
        [sc7] true -> rt: true;
        [sc8] true -> ry: true;
        [sc9] true -> rt: true;
        [sca1] true -> rx: true;
        [sca2] true -> rt: true;
        [sca3] true -> rt: true;



        endmodule""")

    f.close()


"""
Cell Cycle (CC) switch model taken from Cardelli paper
    N: total population size
    ratex: rate for opinion x
    ratey: rate for opinion y
    ratet: rate for all other reactions

    Output: RML file model
"""
def write_cc_model(N, ratex = 1.0, ratey = 1.0, ratet = 1.0):
    # compute remaining number of pure agents
    X = int(1/3 * N)
    f = open(model, "w")
    f.write("""ctmc

        const int N = """ + str(N) + """;
        const int Nhalf = """ + str(X) + """;

        const int xMax = xInit+yInit+bInit;
        const int yMax = xInit+yInit+bInit;
        const int bMax = xInit+yInit+bInit;
        const int wMax = wInit+zInit+uInit;
        const int zMax = wInit+zInit+uInit;
        const int uMax = wInit+zInit+uInit;
        const int pMax = pInit+rInit+qInit;
        const int rMax = pInit+rInit+qInit;
        const int qMax = pInit+rInit+qInit;
        const int sMax = wInit+zInit+uInit;
        const int tMax = pInit+rInit+qInit;

        const int xInit = Nhalf;
        const int yInit = Nhalf;
        const int bInit = Nhalf;
        const int wInit = Nhalf;
        const int zInit = Nhalf;
        const int uInit = Nhalf;
        const int pInit = Nhalf;
        const int rInit = Nhalf;
        const int qInit = Nhalf;
        const int sInit = Nhalf;
        const int tInit = Nhalf;

        module SC

            x : [0..xMax] init xInit;
            y : [0..yMax] init yInit;
            b : [0..bMax] init bInit;
            w : [0..wMax] init wInit;
            z : [0..zMax] init zInit;
            u : [0..uMax] init uInit;
            p : [0..pMax] init pInit;
            r : [0..rMax] init rInit;
            q : [0..qMax] init qInit;
            s : [0..sMax] init sInit;
            t : [0..tMax] init tInit;

            // X + Z -> Z + B
            [cc1] x>0 & z>0 -> (x*z) : (x'=x-1) & (b'=min(b+1, bMax));
            // B + Z -> Z + Y
            [cc2] b>0 & z>0 -> (b*z) : (b'=b-1) & (y'=min(y+1, yMax));
            // Y + R -> R + B
            [cc3] y>0 & r>0 -> (y*r) : (y'=y-1) & (b'=min(b+1, bMax));
            // B + R -> R + X
            [cc4] b>0 & r>0 -> (b*r) : (b'=b-1) & (x'=min(x+1, xMax));

            // W + S -> S + U
            [cc5] w>0 & s>0 -> (w*s) : (w'=w-1) & (u'=min(u+1, uMax));
            // U + S -> S + Z
            [cc6] u>0 & s>0 -> (u*s) : (u'=u-1) & (z'=min(z+1, zMax));
            // Z + X -> X + U
            [cc7] z>0 & x>0 -> (z*x) : (z'=z-1) & (u'=min(u+1, uMax));
            // U + X -> X + W
            [cc8] u>0 & x>0 -> (u*x) : (u'=u-1) & (w'=min(w+1, wMax));
            
            // P + X -> X + Q
            [cc9] p>0 & x>0 -> (p*x) : (p'=p-1) & (q'=min(q+1, qMax));
            // Q + X -> X + R
            [cca1] q>0 & x>0 -> (q*x) : (q'=q-1) & (r'=min(r+1, rMax));
            // R + T -> T + Q
            [cca2] r>0 & t>0 -> (r*t) : (r'=r-1) & (q'=min(q+1, qMax));
            // Q + T -> T + P
            [cca3] q>0 & t>0 -> (q*t) : (q'=q-1) & (p'=min(p+1, pMax));

            [dead1] r=0 & q=0 & x=0 & w=0 & u=0 & b=0 -> (x'=x);

        endmodule

        // base rates
        const double rx = """ + str(ratex/N) + """;
        const double ry = """ + str(ratey/N) + """;
        const double rt = """ + str(ratet/N) + """;

        // module representing the base rates of reactions
        module base_rates
            
            [cc1] true -> ry: true;
            [cc2] true -> ry: true;
            [cc3] true -> rx: true;
            [cc4] true -> rx: true;
            [cc5] true -> rt: true;
            [cc6] true -> ry: true;
            [cc7] true -> rt: true;
            [cc8] true -> rt: true;
            [cc9] true -> rt: true;
            [cca1] true -> rx: true;
            [cca2] true -> rt: true;
            [cca3] true -> rt: true;



        endmodule""")

    f.close()
