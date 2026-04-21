from qfs import qf_isogs,qf_isog_cycle,qf_isog_cycle_power,qf_sibs,cycs_from_ancestors,class_group_id
from nt import discfac,quad_rec,divisors
from identities import *

###########################
# Obtaining generating ls #
###########################

# This is start of bijection algorithm #
def generate_qfs_from_ls(qf0,lset):
    qfs = [qf0]
    for l in lset:
        qfls = []
        for qf in qfs:
            qfls+=qf_isog_cycle(qf,l)
        qfs = qfls
    return qfs

def generate_qfs_from_lks(qf0,lkset):
    qfs = [qf0]
    for lk in lkset:
        qfls = []
        for qf in qfs:
            qfls+=qf_isog_cycle_power(qf,lk)
        qfs = qfls
    return qfs


def qf_cyc_data(d,lcands):
    d0, c = discfac(d)
    qf0 = class_group_id(d)
    lcands = [l for l in lcands if c % l !=0 and quad_rec(d0,l)>=0]
    cycdata = {}
    gens = []
    for l in lcands:
        cycl = qf_isog_cycle(qf0,l)
        if len(cycl)>1 and cycl[1] not in gens and cycl[-1] not in gens:
            cycdata[l] = len(cycl)
            gens.append(cycl[1])
    return cycdata


def qf_cyc_data_ext(d,lcands):
    cycdata = qf_cyc_data(d,lcands)
    ls = [l for l in cycdata]
    for l in ls:
        for k in divisors(cycdata[l])[:-1]:
            cycdata[(l,k)] = cycdata[l]//k
    anc_data= cycs_from_ancestors(class_group_id(d))
    for l in anc_data:
        cycdata[(l,-1)]=len(anc_data[l])
    return {l:cycdata[l] for l in cycdata if type(l)==tuple}



def qf_lk_rig_exts(d,lks,lcands):
    neighbordata = {}
    ls = [lk[0] for lk in lks]
    for l in lcands:
        if quad_rec(d,l) == 1 and l not in ls:
            cyc = qf_isog_cycle(class_group_id(d),l)
            if len(cyc)>1:
                for i in [1,-1]:
                    neighbordata[cyc[i]]=l
    chains = [[class_group_id(d)]]
    for lk in lks:
        cyc = qf_isog_cycle_power(class_group_id(d),lk)
        chains = [ch+[cyc[k]] for ch in chains for k in [1,-1]]
    return {neighbordata[chain[-1]]:chain[-1] for chain in chains if chain[-1] in neighbordata}

def qf_rig_ls_from_ns(d,lcands,ntup):
    # First we figure out size of cycles in each graph
    cycdata = qf_cyc_data_ext(d,lcands)
    n2lks = {n:[] for n in ntup}
    for lk in cycdata:
        if cycdata[lk] in n2lks:
            n2lks[cycdata[lk]].append(lk)
    # We construct tuples of l's that could potentially generate
    ts = [[]]
    cld = 1
    for ni in ntup:
        cld*=ni
        ts = [t+[lk] for t in ts for lk in n2lks[ni] if lk not in t]
    # We now look for a generating set
    gensfd = []
    for ls in ts:
        # Check if it generates
        if len(set(generate_qfs_from_lks(class_group_id(d),ls)))==cld:
            gensfd.append(ls)
            #check if rigid:
            rigexts = qf_lk_rig_exts(d,ls,lcands)
            if len(rigexts)>0:
                return ls, rigexts
    if len(gensfd)>0:
        return gensfd[0],{}
    else:
        return [],{}
            
def qf_ltups_from_ntup(d,lcands,ntup):
    cycdata = qf_cyc_data_ext(d,lcands)
    n2lks = {n:[] for n in ntup}
    for lk in cycdata:
        if cycdata[lk] in n2lks:
            n2lks[cycdata[lk]].append(lk)
    ts = [[]]
    cld = 1
    for ni in ntup:
        cld*=ni
        ts = [t+[lk] for t in ts for lk in n2lks[ni] if lk not in t]
    return [t for t in ts if len(set(generate_qfs_from_lks(class_group_id(d),t)))==cld]

def qf_n2k_check(d,lcands):
    cld = clgr_size_gen(d)
    if abs(d)<13 or cld == 1:
        return [1]
    cycdata = qf_cyc_data_ext(d,lcands)
    # Make sure we got some data
    if len(cycdata)==0:
        raise ValueError('No isogeny data from these primes')
    nmx = max(cycdata.values())
    lks_mx = [lk for lk in cycdata if cycdata[lk]==nmx]
    #If the group is cyclic, we can stop
    if cld == nmx:
        return [lks_mx[0]]
    #We want to check if the group has the form (n,2,2,2,2)
    ntup = [nmx]
    if cld % nmx != 0:
        raise ValueError('Class group order is not divide maximal cycle length, something is wrong') 
    mx_ind = cld//nmx
    while mx_ind % 2 == 0:
        ntup.append(2)
        mx_ind = mx_ind//2
    #When the loop terminates, there are no more potential factors of 2
    #If there's still an odd term, there is no presentation of the desired form
    #so we return the empty list
    if mx_ind>1:
        return []
    n_to_lks = {n:[] for n in ntup}
    for lk in cycdata:
        if cycdata[lk] in n_to_lks:
            n_to_lks[cycdata[lk]].append(lk)
    # We're going to construct all tuples of lks where the first lk has length n
    # and the remaining have length 2
    ts = [[]]
    for ni in ntup:
        ts = [t+[lk] for t in ts for lk in n_to_lks[ni] if lk not in t]
    # We check if any of them is a generating set
    for ls in ts:
        if len(set(generate_qfs_from_lks(class_group_id(d),ls)))==cld:
            return {'ls':ls,'ns':ntup}
    else:
        return {}
    



def qf_riglk_search(d,ntup,lcands):
    cycdata = qf_cyc_data_ext(d,lcands)
    n2lks = {n:[] for n in ntup}
    for lk in cycdata:
        if cycdata[lk] in n2lks:
            n2lks[cycdata[lk]].append(lk)
    ts = [[]]
    cld = 1
    for ni in ntup:
        cld*=ni
        ts = [t+[lk] for t in ts for lk in n2lks[ni] if lk not in t]
    gensfound = []
    for lks in ts:
        # Check if they generate
        if len(set(generate_qfs_from_lks(class_group_id(d),lks)))==cld:
            # They do generate; check orientation
            rigends = qf_lk_rig_exts(d,lks,lcands)
            if len(rigends)>0:
                return {'lks':lks,'ns':ntup,'rigext':rigends}
            else:
                gensfound.append(lks)
    if len(gensfound)>0:
        return {'lks':gensfound[0],'ns':ntup,'rigext':{}}
    else:
        return {'lks':[],'ns':ntup,'rigext':{}}
    
def qf_lgen_search(d,lcands,md=None,checkn2 = True):
    cld = clgr_size_gen(d)
    if abs(d)<13 or cld==1:
        return [1]
    if checkn2:
        n2ch = qf_n2k_check(d,lcands)
        if len(n2ch)>0:
            return n2ch
    cycdata = qf_cyc_data_ext(d,lcands)
    if md == None:
        md = len(cycdata)
    # We don't have a cyclic group, so we will need multiple generators
    # Before looking for a generating set we determine whether we will also
    # need an orientation
    lks_scores = {(lk,):cycdata[lk] for lk in cycdata}
    for i in range(max(md,len(lcands))):
        newdata = {}
        for lks in lks_scores:
            ls = [lk[0] for lk in lks]
            for lk1 in cycdata:
                lks1_score = (lks_scores[lks]*cycdata[lk1])
                if lk1[0] > max(ls) and cld % lks1_score==0:
                    lks1_l = list(lks)+[lk1]
                    actual_score = len(set(generate_qfs_from_lks(class_group_id(d),lks1_l)))
                    if lks1_score == actual_score:
                        newdata[tuple(lks1_l)]=lks1_score
                        if lks1_score == cld:
                            ntup = tuple([cycdata[lk] for lk in lks1_l])
                            # We have a generating set! 
                            # Now just check orientation
                            rigdata = qf_lk_rig_exts(d,lks1_l,lcands)
                            if len(rigdata)>0:
                                return {'ls':lks1_l,'ns':ntup,'rigext':rigdata}
                            else:
                                return qf_riglk_search(d,ntup,lcands)
        if len(newdata)>0:
            lks_scores = newdata
        else:
            ms = max(lks_scores.values())
            lks_mx = [lk for lk in lks_scores if lks_scores[lk]==ms]
            return {'ls':tuple(lks_mx[0]),'ns':(),'rigext':{}}
    ms = max(lks_scores.values())
    lks_mx = [lk for lk in lks_scores if lks_scores[lk]==ms]
    return {'ls':tuple(lks_mx[0]),'ns':(),'rigext':{}}



