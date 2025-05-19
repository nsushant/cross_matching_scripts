import pynbody 
import darklight 
import pandas as pd 
import numpy as np 
import tangos 
from tangos.examples.mergers import *
from particle_tagging.edge.utils import *
import sys 
import os
from scipy.spatial.distance import cdist
### Function defs 
pynbody.config["halo-class-priority"] = [pynbody.halo.hop.HOPCatalogue]

def calculate_Rizzo_energy_distance(X,Y):

    # Returns a measure of simmilarity between two distributions called the Energy Distance 
    # A lower energy distance = distributions are more simmilar (vice-versa)
    # citation: WIREs Comput Stat 2016, 8:27–38. doi: 10.1002/wics.1375 (see sec. TESTING FOR EQUAL DISTRIBUTIONS)
    
    XYMeanCdist = np.mean(cdist(X, Y))   # Spread between the two distributions of samples
    XXMeanCdist = np.mean(cdist(X, X))   # Spread within distribution X
    YYMeanCdist = np.mean(cdist(Y, Y))   # Spread within distribution Y
    
    return 2 * XYMeanCdist - XXMeanCdist - YYMeanCdist 
    


def group_mergers(z_merges,h_merges,q_merges):
    #groups the halo objects of merging halos by redshift                                                                                                   
    # quantities the function outputs 
    merging_halos_grouped_by_z = []
    qvalues_grouped_by_z = []
    m200_groupd_by_z = []
    z_unique_values = sorted(list(set(z_merges)))
    
    #internal array used to avoid halo overlaps
    list_of_selected_tuples = []

    # for each of these redshifts                                                                                                                
    for i in z_unique_values:
        # indices of halo objects merging at given redshift 'i'                                                                                              
        # zmerges = 1D (can contain 'i' multiple times)                                                                                          
        lists_of_halos_merging_at_current_z = np.where(z_merges==i)


        all_halos_merging_at_current_z=[]
        qvalsz = []
        m200s_list_main = []

        # collect halo objects of merging halos for the redhsift 'i'                                                                                   
        for list_of_halos in lists_of_halos_merging_at_current_z :
            halos_merging_at_current_z =np.array([])
            qvals = np.array([])
            m200s_list = np.array([])

            for merging_halo_object in list_of_halos:

                #halos_merging_at_current_z = np.append(halos_merging_at_current_z, h_merges[merging_halo_object][1:])                                       
                #print("len hmerge: ",len(h_merges[merging_halo_object][1:]))                                                                                
                #qvals = np.append(qvals,q_merges[merging_halo_object])                                                                                      
                haloM =  h_merges[merging_halo_object][1:][0]

                halonums_over_life, ts_over_life = haloM.calculate_for_progenitors("halo_number()","t()")

                overlap = False


                #print( halonums_over_life, ts_over_life)                                                                                                    
                for h,t in zip(halonums_over_life,ts_over_life):

                    pair = [h,t]

                    if pair in list_of_selected_tuples:
                        print("overlap found")
                        overlap = True
                        

                    if (overlap == False):

                        list_of_selected_tuples.append(pair)

                if (overlap==True):
                    continue 

                print("FALSE OVERLAP")
                halos_merging_at_current_z = np.append(halos_merging_at_current_z, h_merges[merging_halo_object][1:])

                qvals = np.append(qvals,q_merges[merging_halo_object])

                if np.isin("M200c_DM",haloM.keys()):
                    #print(haloM["M200c"])                                                                                                                   
                    m200s_list = np.append(m200s_list,haloM["M200c_DM"])

                else:
                    m200s_list = np.append(m200s_list,0)

                #halonums_over_life, ts_over_life = haloM.calculate_for_progenitors("halo_number()","t()")                                                   

                #for h,t in zip(halonums_over_life[0],ts_over_life[0]):                                                                                      

                #if np.logical_not(np.isin(tuple(h,t),lists_of_halos_merging_at_current_z)):                                                                 

                #lists_of_halos_merging_at_current_z = np.append(lists_of_halos_merging_at_current_z,tuple(h,t))                                             


            all_halos_merging_at_current_z.append(halos_merging_at_current_z)

            qvalsz.append(qvals)

            m200s_list_main.append(m200s_list)


        merging_halos_grouped_by_z.append(all_halos_merging_at_current_z)
        qvalues_grouped_by_z.append(qvalsz)

        m200_groupd_by_z.append(m200s_list_main)

    return merging_halos_grouped_by_z, z_unique_values, qvalues_grouped_by_z, m200_groupd_by_z


def EuclideanDistance(xyz1,xyz2):

    x = xyz1[0] - xyz2[0]
    y = xyz1[1] - xyz2[1]
    z = xyz1[2] - xyz2[2]

    return np.sqrt( x**2 + y**2 + z**2 )


### Input Processing

haloname = str(sys.argv[1])

if len(str(haloname)) <= 8:
    DMOname = haloname+"_DMO"
    HYDROname = haloname+"_fiducial"

else: 
    HaloNameDecomp = haloname.split("_")
    DMOname = HaloNameDecomp[0]+"_DMO_"+HaloNameDecomp[2]
    HYDROname = HaloNameDecomp[0]+"_fiducial_"+HaloNameDecomp[2]
    haloname = HaloNameDecomp[0]


pynbody_path = '/scratch/dp101/shared/EDGE/'
# Finding best match in tangos db

# Start of crossreff 
tangos.core.init_db("/scratch/dp101/shared/EDGE/tangos/"+str(haloname)+".db")

## Get DMO data 
DMOsim = tangos.get_simulation(DMOname)
DMOMain = DMOsim.timesteps[-1].halos[0]
HaloNumsDMO = DMOMain.calculate_for_progenitors("halo_number()")[0][::-1]
RedsMainDMO = DMOMain.calculate_for_progenitors("z()")[0][::-1]
# create arrays that associate tangos timesteps with stored redshifts 
RedshiftsDMO= []
TimeStepIdxsDMO = []
tstepidx = 0
for t in range(len(DMOsim.timesteps[:])):
    
    if len(DMOsim.timesteps[t].halos[:]) == 0: 
        tstepidx+=1
        continue 
    else: 
        RedshiftsDMO.append(DMOsim.timesteps[t].halos[0].calculate("z()"))
        TimeStepIdxsDMO.append(tstepidx)
        tstepidx+=1

RedshiftsDMO = np.asarray(RedshiftsDMO)
TimeStepIdxsDMO = np.asarray(TimeStepIdxsDMO)

# these two arrays should have the same length
print("DMO:",len(RedshiftsDMO),len(HaloNumsDMO))

## Get HYDRO data 
HYDROsim = tangos.get_simulation(HYDROname)
HYDROMain = HYDROsim.timesteps[-1].halos[0]
HaloNumsHYDRO = HYDROMain.calculate_for_progenitors("halo_number()")[0][::-1]
RedshiftsTangosMainHYDRO = HYDROMain.calculate_for_progenitors("z()")[0][::-1]
RedshiftsHYDRO = []
TimeStepIdxsHYDRO = []
TimeStepsHYDRO = []
tstepidx = 0
for t in range(len(HYDROsim.timesteps[:])):
    
    if (len(HYDROsim.timesteps[t].halos[:]) == 0):
        tstepidx+=1
        continue
    else:
        RedshiftsHYDRO.append(HYDROsim.timesteps[t].halos[0].calculate("z()"))
        TimeStepIdxsHYDRO.append(tstepidx)
        TimeStepsHYDRO.append(str(HYDROsim.timesteps[t]))
        tstepidx+=1


TimeStepsHYDRO = np.asarray(TimeStepsHYDRO)
TimeStepIdxsHYDRO = np.asarray(TimeStepIdxsHYDRO)
RedshiftsHYDROAll = np.asarray(RedshiftsHYDRO)

# these should have the same lengths
print("HYDRO:",len(HaloNumsHYDRO),len(RedshiftsHYDRO))

#Processing Merger Tree 
MergerRedshiftsHYDRO, MergerRatiosHYDRO, MergerHaloObjectsHYDRO = get_mergers_of_major_progenitor(HYDROMain)
GroupedHalosHYDRO, GroupedRedshiftsHYDRO, GroupedMergerRatiosHYDRO,Groupedm200s = group_mergers(MergerRedshiftsHYDRO, MergerHaloObjectsHYDRO, MergerRatiosHYDRO)

idx_of_best_match_DMO = [np.argmin(RedshiftsDMO[np.where(RedshiftsDMO > zh)]) for zh in GroupedRedshiftsHYDRO]

#tstepidxsHYDRO = TimeStepIdxsHYDRO[np.asarray(idx_of_best_match_hydro)] 
tstepidxsDMO = TimeStepIdxsDMO[np.asarray(idx_of_best_match_DMO)]


#print(idx_of_best_match_hydro)
hydrohalo_matched = []
dmohalo_matched = [] 
HydroHaloMstars = []

for z in range(len(GroupedRedshiftsHYDRO))[::-1]:

    HYDROMergingHalosThisRedshift = GroupedHalosHYDRO[z][0]            
    if (len(HYDROMergingHalosThisRedshift) == 0): 
        continue

    MergerTimestep = HYDROMergingHalosThisRedshift[0].timestep

    HYDROTimestepThisMerger = np.where(TimeStepsHYDRO == str(MergerTimestep))[0][0]
    
    HYDROhalonumidx = np.where(RedshiftsTangosMainHYDRO == RedshiftsHYDRO[HYDROTimestepThisMerger])[0][0]
    
    HYDROMainHaloThisRedshift = HYDROsim.timesteps[ TimeStepIdxsHYDRO[HYDROTimestepThisMerger] ].halos[ int(HaloNumsHYDRO[HYDROhalonumidx]) - 1 ]

    DMOhalonumidx = np.where(RedsMainDMO==RedshiftsDMO[tstepidxsDMO[z]])[0][0]
    MainHaloDMOThisRedshift = DMOsim.timesteps[ tstepidxsDMO[z] ].halos[ int(HaloNumsDMO[DMOhalonumidx]) - 1 ]
    DMOHalosThisRedshift = list(DMOsim.timesteps[ tstepidxsDMO[z] ].halos[:])

    DMOHalosThisRedshift.remove(MainHaloDMOThisRedshift)
    dm_mass = []    

    for DMOhalo in DMOHalosThisRedshift:
            
        try: 
            dm_mass.append(MainHaloDMOThisRedshift.calculate("M200c")/DMOhalo.calculate("M200c"))

        except:
            #print(e)
            dm_mass.append(0)

    # load in HYDRO data 
    outputHYDRO = str(HYDROsim.timesteps[ TimeStepIdxsHYDRO[HYDROTimestepThisMerger]  ]).split("/")[-1]
    simfnHYDRO = os.path.join(pynbody_path,HYDROname,outputHYDRO)
    HYDROParticles = pynbody.load(simfnHYDRO)
    pynbody.analysis.halo.center(HYDROParticles.halos(int(HaloNumsHYDRO[HYDROhalonumidx]) - 1))
    HYDROParticles.physical_units()
    
    ParticlesLoadedIn = False
    
    for MergingHYDROhalo in HYDROMergingHalosThisRedshift:
        print(MergingHYDROhalo)
        try: 
            MergingHYDROhalo["M200c_stars"]
            
            if MergingHYDROhalo["M200c_stars"] == 0: 
                print("No Mstar")
                continue 
            
        except Exception as e: 
            print(e)
            continue

        
        try:
            m200MergingHYDROhalo = HYDROMainHaloThisRedshift["M200c"]/MergingHYDROhalo["M200c"]
        
            print("added")

        except Exception as er:
            print(er)
            continue


        # sorts mass difference in M200 in ascending order    
        closest_mass_match = np.argsort(np.abs(np.asarray(dm_mass) - m200MergingHYDROhalo))[:5]                                  

        ### Tangos part ends and energy distance calculation starts 
        #HYDRO halo 6D array
        MergingHYDROHaloNumber = MergingHYDROhalo.calculate("halo_number()")
        MergingHYDROHaloParticles = HYDROParticles.halos()[int(MergingHYDROHaloNumber)-1]

        px = MergingHYDROHaloParticles.dm["vel"][:,0]*MergingHYDROHaloParticles.dm["mass"]
        py = MergingHYDROHaloParticles.dm["vel"][:,1]*MergingHYDROHaloParticles.dm["mass"]
        pz = MergingHYDROHaloParticles.dm["vel"][:,2]*MergingHYDROHaloParticles.dm["mass"]
        
        PhaseArrayHydro = np.stack((MergingHYDROHaloParticles.d['x'], MergingHYDROHaloParticles.d['y'], MergingHYDROHaloParticles.d['z'], px, py, pz), axis=1)

        # Load in DMO data
        if ParticlesLoadedIn == False:          
            # load in DMO pynbody data
            output_num = str(DMOsim.timesteps[ tstepidxsDMO[z] ]).split("/")[-1]
            simfnDMO = os.path.join(pynbody_path,DMOname,output_num)
            DMOParticles = pynbody.load(simfnDMO)


        PhaseArraysDMO = []
        for MassMatch in closest_mass_match:
            MassMatchedDMOHalo = DMOHalosThisRedshift[MassMatch]
            HaloNumMassMatchedDMOHalo = MassMatchedDMOHalo.calculate("halo_number()")
            MassMatchedDMOHaloParticles = DMOParticles.halos()[int(HaloNumMassMatchedDMOHalo) - 1]

            pxd = MassMatchedDMOHaloParticles["vel"][:,0]*MassMatchedDMOHaloParticles["mass"]
            pyd = MassMatchedDMOHaloParticles["vel"][:,1]*MassMatchedDMOHaloParticles["mass"]
            pzd = MassMatchedDMOHaloParticles["vel"][:,2]*MassMatchedDMOHaloParticles["mass"]
            PhaseArraysDMO.append(np.stack((MassMatchedDMOHaloParticles["x"],MassMatchedDMOHaloParticles["y"],MassMatchedDMOHaloParticles["z"],pxd,pyd,pzd),axis=1))
            

        EnergyDistances = np.array([])

        for PhaseArrayDMO in PhaseArraysDMO: 

            EnergyDistanceHalo = calculate_Rizzo_energy_distance(PhaseArrayDMO,PhaseArrayHydro)
            EnergyDistances = np.append(EnergyDistances,EnergyDistanceHalo)

        
        best_match_2_fold =  closest_mass_match[np.argmin(EnergyDistances)]
        
        hydrohalo_matched.append(MergingHYDROhalo)
        HydroHaloMstars.append(MergingHYDROhalo["M200c_stars"])
        dmohalo_matched.append(DMOHalosThisRedshift[best_match_2_fold])
        #dmohalo_matched.append(d.timesteps[idx_of_best_match[z]].halos[closest_mass_match[0]])

print(hydrohalo_matched)
print(dmohalo_matched)

df = pd.DataFrame({"halo":dmohalo_matched,"mstar":HydroHaloMstars,"hydrohalo":hydrohalo_matched})                                                                  

df.to_csv("dmo_hydro_crossreffs/TwoFoldCrossreff_"+DMOname+".csv")                    
