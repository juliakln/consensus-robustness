"""
Define properties in BLTL for Plasmalab.
Reaching a stable consensus, switching consensus, recovery analysis. 

Save in '../models/consensus.bltl"

"""

import os


# change to location of this script
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

property = "../models/consensus.bltl"



"""
Write property for reaching stable consensus with stubborn individuals for Plasmalab
    N: total population size
    stubborn: z (zealots), c (contrarians) or b (both)
    majority: more than m% of population commits to same decision
    distance: difference of at least d between majority and those favouring opposite decision
    reaching: consensus is achieved within t minutes from start of dynamics
    holding: group maintains consensus for at least h minutes

    Output: BLTL file property
"""
def write_property_stableconsensus(N = 100, stubborn = 'z', majority = 50, distance = 10, reaching = 35, holding = 40):
    # compute absolute number to reach majority
    threshold = int((majority / 100) * N)
    f = open(property, "w")
    if stubborn == 'z':
        f.write("F<="+str(reaching)+" (G<="+str(holding)+" (((x+Zx>="+str(threshold)+") & (x-y>="+str(distance)+")) | ((y+Zy>="+str(threshold) + ") & (y-x>="+str(distance)+"))))")
    elif stubborn == 'c':
        f.write("F<="+str(reaching)+" (G<="+str(holding)+" (((x+Cx>="+str(threshold)+") & ((x+Cx)-(y+Cy)>="+str(distance)+")) | ((y+Cy>="+str(threshold) + ") & ((y+Cy)-(x+Cx)>="+str(distance)+"))))")
    elif stubborn == 'b':
        f.write("F<="+str(reaching)+" (G<="+str(holding)+" (((x+Cx+Zx>="+str(threshold)+") & ((x+Cx+Zx)-(y+Cy+Zy)>="+str(distance)+")) | ((y+Cy+Zy>="+str(threshold) + ") & ((y+Cy+Zy)-(x+Cx+Zx)>="+str(distance)+"))))")
    else:
        raise Exception('Type of stubborn individual not supported.')
    f.close()


"""
Reaching a stable consensus for option X only
    N: total population size
    stubborn: z (zealots), c (contrarians), b (both) or none
    majority: more than m% of population commits to same decision
    distance: difference of at least d between majority and those favouring opposite decision
    reaching: consensus is achieved within t minutes from start of dynamics
    holding: group maintains consensus for at least h minutes
    
    Output: BLTL file property
"""
def write_property_stableconsensus_asym_x(N = 100, stubborn = 'z', majority = 50, distance = 10, reaching = 35, holding = 40):
    # compute absolute number to reach majority
    threshold = int((majority / 100) * N)
    f = open(property, "w")
    if stubborn == 'z':
        f.write("F<="+str(reaching)+" (G<="+str(holding)+" (((x+Zx>="+str(threshold)+") & (x-y>="+str(distance)+"))))")
    elif stubborn == 'c':
        f.write("F<="+str(reaching)+" (G<="+str(holding)+" (((x+Cx>="+str(threshold)+") & ((x+Cx)-(y+Cy)>="+str(distance)+"))))")
    elif stubborn == 'b':
        f.write("F<="+str(reaching)+" (G<="+str(holding)+" (((x+Cx+Zx>="+str(threshold)+") & ((x+Cx+Zx)-(y+Cy+Zy)>="+str(distance)+"))))")
    elif stubborn == 'none':
        f.write("F<="+str(reaching)+" (G<="+str(holding)+" (((x>="+str(threshold)+") & ((x-y)>="+str(distance)+"))))")
    else:
        raise Exception('Type of stubborn individual not supported.')
    f.close()

"""
Reaching a stable consensus for option Y only
"""
def write_property_stableconsensus_asym_y(N = 100, stubborn = 'z', majority = 50, distance = 10, reaching = 35, holding = 40):
    # compute absolute number to reach majority
    threshold = int((majority / 100) * N)
    f = open(property, "w")
    if stubborn == 'z':
        f.write("F<="+str(reaching)+" (G<="+str(holding)+" (((y+Zy>="+str(threshold)+") & ((y+Zy)-(x+Zx)>="+str(distance)+"))))")
    elif stubborn == 'c':
        f.write("F<="+str(reaching)+" (G<="+str(holding)+" (((y+Cy>="+str(threshold)+") & ((y+Cy)-(x+Cx)>="+str(distance)+"))))")
    elif stubborn == 'b':
        f.write("F<="+str(reaching)+" (G<="+str(holding)+" (((y+Cy+Zy>="+str(threshold)+") & ((y+Cy+Zy)-(x+Cx+Zx)>="+str(distance)+"))))")
    else:
        raise Exception('Type of stubborn individual not supported.')
    f.close()


"""
Write property for switching consensus 
    stubborn: z (zealots), c (contrarians) or b (both)
    distance: difference of at least d between majority and those favouring opposite decision
    reaching: consensus is achieved within t minutes from start of dynamics
    holding: group maintains consensus for at least h minutes

    Output: BLTL file property
"""
def write_property_switching(N, stubborn = 'z', majority = 50, distance = 10, reaching = 35, holding = 10):
    f = open(property, "w")
    if stubborn == 'z':
        f.write("F<="+str(reaching)+" (((x+Zx)-(y+Zy)>="+str(distance)+" & (true U<="+str(holding)+" ((y+Zy)-(x+Zx)>="+str(distance)+"))) | ((y+Zy)-(x+Zx)>="+str(distance)+" & (true U<="+str(holding)+" ((x+Zx)-(y+Zy)>="+str(distance)+"))))")
    elif stubborn == 'c':
        f.write("F<="+str(reaching)+" (((x+Cx)-(y+Cy)>="+str(distance)+" & (true U<="+str(holding)+" ((y+Cy)-(x+Cx)>="+str(distance)+"))) | ((y+Cy)-(x+Cx)>="+str(distance)+" & (true U<="+str(holding)+" ((x+Cx)-(y+Cy)>="+str(distance)+"))))")
    elif stubborn == 'b':
        f.write("F<="+str(reaching)+" (((x+Cx)-(y+Cy)>="+str(distance)+" & (true U<="+str(holding)+" ((y+Cy)-(x+Cx)>="+str(distance)+"))) | ((y+Cy)-(x+Cx)>="+str(distance)+" & (true U<="+str(holding)+" ((x+Cx)-(y+Cy)>="+str(distance)+"))))")
    else:
        raise Exception('Type of stubborn individual not supported.')
    f.close()


"""
Write property for recovery for Plasmalab
    N: total population size
    stubborn: NOTE so far only no disruptives supported
    majority: more than m% of population commits to same decision
    distance: difference of at least d between majority and those favouring opposite decision
    G: holding time, group maintains consensus for at least h minutes
    rec_time: time to recover from majority for worse option

    Output: BLTL file property
"""
def write_property_recovery(N = 100, stubborn = 'z', majority = 50, distance = 10, reaching = 10, G = 100, rec_time = 10):
    threshold = int((majority / 100) * N)
    f = open(property, "w")
    if stubborn == 'z':
        f.write("(! ((y+Zy>="+str(threshold)+") & ( (y+Zy)-(x+Zx)>="+str(distance)+") )) U<="+str(G)+" ( ((y+Zy>="+str(threshold)+") & ( (y+Zy)-(x+Zx)>="+str(distance)+")) & (F<="+str(rec_time)+" ( G<=100 ( ! (( y+Zy>="+str(threshold)+") & ( (y+Zy)-(x+Zx)>="+str(distance)+")) ) ) ) )")
    else:
        raise Exception('Type of stubborn individual not supported.')
    f.close()


"""
Recovery from Y to good majority for X
"""
def write_property_recoverygood(N = 100, stubborn = 'z', majority = 50, distance = 10, F = 10, G = 100, rec_time = 10):
    threshold = int((majority / 100) * N)
    f = open(property, "w")
    if stubborn == 'z':
        f.write("(! ((y>="+str(threshold)+") & ( (y-x)>="+str(distance)+") )) U<="+str(F)+" ( ((y>="+str(threshold)+") & ( (y-x)>="+str(distance)+")) & (F<="+str(rec_time)+" ( G<=" + str(G) + " ( (( x>="+str(threshold)+") & ( (x-y)>="+str(distance)+")) ) ) ) )")
    elif stubborn == 'c':
        f.write("G<="+str(G)+" ( !( (y+Cy>="+str(threshold)+") & ( (y+Cy)-(x+Cx)>="+str(distance)+") ) | ( F<="+str(rec_time)+" ( G<=100 ( !( (y+Cy>="+str(threshold)+") & ( (y+Cy)-(x+Cx)>="+str(distance)+") ) ) ) ) )")
    else:
        raise Exception('Type of stubborn individual not supported.')
    f.close()


"""
Reach consensus for X and recover from majority for Y
"""
def write_property_recoveryandconsensus(N = 100, stubborn = 'z', majority = 50, distance = 10, F = 10, G = 50, rec_time = 1):
    threshold = int((majority / 100) * N)
    f = open(property, "w")
    if stubborn == 'z':
        f.write("(F<="+str(F)+" (G<="+str(G)+" (((x>="+str(threshold)+") & ((x-y)>="+str(distance)+"))))) | ((! ((y+Zy>="+str(threshold)+") & ( (y+Zy)-(x+Zx)>="+str(distance)+") )) U<="+str(F)+" ( ((y+Zy>="+str(threshold)+") & ( (y+Zy)-(x+Zx)>="+str(distance)+")) & (F<="+str(rec_time)+" ( G<="+str(G)+" ( ! (( y+Zy>="+str(threshold)+") & ( (y+Zy)-(x+Zx)>="+str(distance)+")) ) ) ) ))")
    else:
        raise Exception('Type of stubborn individual not supported.')
    f.close()


"""
Reach consensus for X and recover from majority for Y to majority for X
"""
def write_property_recoveryandconsensusgood(N = 100, stubborn = 'z', majority = 50, distance = 10, F = 10, G = 50, rec_time = 1):
    threshold = int((majority / 100) * N)
    f = open(property, "w")
    if stubborn == 'z':
        f.write("(F<="+str(F)+" (G<="+str(G)+" (((x>="+str(threshold)+") & ((x-y)>="+str(distance)+"))))) | ((! ((y>="+str(threshold)+") & ( (y-x)>="+str(distance)+") )) U<="+str(F)+" ( ((y>="+str(threshold)+") & ( (y-x)>="+str(distance)+")) & (F<="+str(rec_time)+" ( G<="+str(G)+" ( (( x>="+str(threshold)+") & ( (x-y)>="+str(distance)+")) ) ) ) ))")
    else:
        raise Exception('Type of stubborn individual not supported.')
    f.close()

