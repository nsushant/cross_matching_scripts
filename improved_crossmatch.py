import pynbody 
import darklight 
import pandas as pd 
import numpy as np 
import tangos 
from tangos.examples.mergers import *

#from cross_matching_hydro_and_dmo import group_mergers
import sys 



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


haloname = str(sys.argv[1])

# Getting inputs
if len(str(haloname)) <= 8:
    DMOname = haloname+"_DMO"
    HYDROname = haloname+"_fiducial"

else: 
    HaloNameDecomp = haloname.split("_")
    DMOname = HaloNameDecomp[0]+"_DMO_"+HaloNameDecomp[2]
    HYDROname = HaloNameDecomp[0]+"_fiducial_"+HaloNameDecomp[2]
    haloname = HaloNameDecomp[0]

tangos.core.init_db("/scratch/dp101/shared/EDGE/tangos/"+str(haloname)+".db")
# DMO data 
DMOsim = tangos.get_simulation(DMOname)
DMOMain = DMOsim.timesteps[-1].halos[0]
HaloNumsDMO = DMOMain.calculate_for_progenitors("halo_number()")[0][::-1]
RedshiftsDMOAll= []
TimeStepIdxsDMO = []
print(DMOsim.timesteps[0].halos[:])
tstepidx = 0
for t in range(len(DMOsim.timesteps[:])):
    
    if len(DMOsim.timesteps[t].halos[:]) == 0: 
        RedshiftsDMOAll.append(99999)
        TimeStepIdxsDMO.append(tstepidx)
        tstepidx+=1
        continue 
    else: 
        RedshiftsDMOAll.append(DMOsim.timesteps[t].halos[0].calculate("z()"))
        TimeStepIdxsDMO.append(tstepidx)
        tstepidx+=1

TimeStepIdxsDMO = np.asarray(TimeStepIdxsDMO)
RedshiftsDMOAll = np.asarray(RedshiftsDMOAll)
RedshiftsTangosMainDMO = DMOMain.calculate_for_progenitors("z()")[0][::-1]
#print(np.isin(RedshiftsDMOAll,RedshiftsTangosMainDMO))
RedshiftsDMO = RedshiftsDMOAll[np.isin(RedshiftsDMOAll,np.asarray(RedshiftsTangosMainDMO))]
TimeStepIdxsDMO = TimeStepIdxsDMO[np.isin(RedshiftsDMOAll,np.asarray(RedshiftsTangosMainDMO))]

print("DMO:",len(RedshiftsDMO),len(HaloNumsDMO))



# HYDRO data                                                                                                                         
HYDROsim = tangos.get_simulation(HYDROname)
HYDROMain = HYDROsim.timesteps[-1].halos[0]
HaloNumsHYDRO = HYDROMain.calculate_for_progenitors("halo_number()")[0][::-1]
RedshiftsTangosMainHYDRO = HYDROMain.calculate_for_progenitors("z()")[0][::-1]
RedshiftsHYDROAll = []
TimeStepIdxsHYDRO = []
tstepidx = 0
for t in range(len(HYDROsim.timesteps[:])):
    
    if len(HYDROsim.timesteps[t].halos[:]) == 0:
        RedshiftsHYDROAll.append(99999)
        TimeStepIdxsHYDRO.append(tstepidx)
        tstepidx+=1
        continue
    else:
        RedshiftsHYDROAll.append(HYDROsim.timesteps[t].halos[0].calculate("z()"))
        TimeStepIdxsHYDRO.append(tstepidx)
        tstepidx+=1



TimeStepIdxsHYDRO = np.asarray(TimeStepIdxsHYDRO)
RedshiftsHYDROAll = np.asarray(RedshiftsHYDROAll)
RedshiftsHYDRO = RedshiftsHYDROAll[np.isin(RedshiftsHYDROAll,np.asarray(RedshiftsTangosMainHYDRO))]
TimeStepIdxsHYDRO = TimeStepIdxsHYDRO[np.isin(RedshiftsHYDROAll,np.asarray(RedshiftsTangosMainHYDRO))]
print("HYDRO:",len(HaloNumsHYDRO),len(RedshiftsHYDRO))

#Processing Merger Tree 
MergerRedshiftsHYDRO, MergerRatiosHYDRO, MergerHaloObjectsHYDRO = get_mergers_of_major_progenitor(HYDROMain)
GroupedHalosHYDRO, GroupedRedshiftsHYDRO, GroupedMergerRatiosHYDRO,Groupedm200s = group_mergers(MergerRedshiftsHYDRO, MergerHaloObjectsHYDRO, MergerRatiosHYDRO)
simpath = "/scratch/dp101/shared/EDGE/"

print(RedshiftsHYDRO,HaloNumsHYDRO)

#idx_of_best_match_DMOHnums = [np.argmin(abs(RedshiftsDMO - zh)) for zh in GroupedRedshiftsHYDRO]
#idx_of_best_match_hydroHnums = [np.argmin(abs(RedshiftsHYDRO - zh)) for zh in GroupedRedshiftsHYDRO]



idx_of_best_match_DMO = [np.argmin(abs(RedshiftsDMO - zh)) for zh in GroupedRedshiftsHYDRO]
idx_of_best_match_hydro = [np.argmin(abs(RedshiftsHYDRO - zh)) for zh in GroupedRedshiftsHYDRO]

tstepidxsHYDRO = TimeStepIdxsHYDRO[np.asarray(idx_of_best_match_hydro)] 
tstepidxsDMO = TimeStepIdxsDMO[np.asarray(idx_of_best_match_DMO)]


print(idx_of_best_match_hydro)

hydrohalo_matched = []
dmohalo_matched = [] 
HydroHaloMstars = []

for z in range(len(GroupedRedshiftsHYDRO))[::-1]:

    #index_z_hydro = np.where(zhydro>=zdmo[z])[-1]                                                                                      
    
    #print("index-->",index_z_hydro)                                                                                                                 
    #print(np.where(zhydro>=zdmo[z]))                                                                                                                
    HYDROMergingHalosThisRedshift = GroupedHalosHYDRO[z][0]
    #print("hydrohalo:",haloshydro.calculate("output()"))                                                                         
    #print(idx_of_best_match_hydro[z],len(HYDROsim.timesteps[ idx_of_best_match_hydro[z] - 1 ].halos[:]))
    MainHaloHYDROThisRedshift = HYDROsim.timesteps[ tstepidxsHYDRO[z] - 3].halos[ int(HaloNumsHYDRO[idx_of_best_match_hydro[z] - 3]) - 1 ]

    #print(MainHaloHYDROThisRedshift,HYDROMergingHalosThisRedshift)
    
    MainHaloDMOThisRedshift = DMOsim.timesteps[ tstepidxsDMO[z] - 1].halos[ int(HaloNumsDMO[ idx_of_best_match_DMO[z] - 1]) - 1 ]
    DMOHalosThisRedshift = list(DMOsim.timesteps[ tstepidxsDMO[z] - 1].halos[:])
    #[int(HaloNumsDMO[ idx_of_best_match_DMO[z]]):]

    DMOHalosThisRedshift.remove(MainHaloDMOThisRedshift)
    dm_mass = []    

    for DMOhalo in DMOHalosThisRedshift:
            
        try: 
            dm_mass.append(DMOhalo.calculate("M200c"))

        except:
            dm_mass.append(0)

            
    
    for MergingHYDROhalo in HYDROMergingHalosThisRedshift:
        print(MergingHYDROhalo)
        try: 
            MergingHYDROhalo["M200c_stars"]
            
            if MergingHYDROhalo["M200c_stars"] == 0: 
                continue 
            
        except: 
            continue


        try:
            #mainhalo and merging halo dist in hydro sim

            DistanceFromMainHYDROHalo = EuclideanDistance(np.asarray(MainHaloHYDROThisRedshift["shrink_center"]),np.asarray(MergingHYDROhalo.calculate("shrink_center")))

            m200MergingHYDROhalo = MergingHYDROhalo["M200c_DM"]
        
            
            print("added")

        except:
            continue


        # sorts mass difference in M200 in ascending order    
        closest_mass_match = np.argsort(np.abs(np.asarray(dm_mass) - m200MergingHYDROhalo))[:10]
        print(closest_mass_match)
        
        #print("closest match:",np.log10(np.abs(dm_mass[closest_mass_match])),np.log10(m200MergingHYDROhalo))

        DistancesFromMainDMOHalo = []
        
        for MassMatch in closest_mass_match:

            try:
                #MassMatchCen = DMOHalosThisRedshift[MassMatch].calculate("shrink_center")
                
                
                MassMatchCen = DMOHalosThisRedshift[MassMatch].calculate("shrink_center")
            
                dist = EuclideanDistance(np.asarray(MainHaloDMOThisRedshift["shrink_center"]),np.asarray(MassMatchCen))
        
                print(DMOHalosThisRedshift[MassMatch],dist,DistanceFromMainHYDROHalo,MainHaloHYDROThisRedshift)
                
                DistancesFromMainDMOHalo.append(dist)
            
            except Exception as e:
                print(e)
                DistancesFromMainDMOHalo.append(999999)
                continue 
            
        best_match_2_fold = np.argmin(np.abs(np.asarray(DistancesFromMainDMOHalo)-DistanceFromMainHYDROHalo))

        hydrohalo_matched.append(MergingHYDROhalo)
        
        HydroHaloMstars.append(MergingHYDROhalo["M200c_stars"])
        dmohalo_matched.append(DMOHalosThisRedshift[closest_mass_match[int(best_match_2_fold)]])
        #dmohalo_matched.append(d.timesteps[idx_of_best_match[z]].halos[closest_mass_match[0]])

print(hydrohalo_matched)
print(dmohalo_matched)

df = pd.DataFrame({"halo":dmohalo_matched,"mstar":HydroHaloMstars,"hydrohalo":hydrohalo_matched})                                                                  


df.to_csv("dmo_hydro_crossreffs/TwoFoldCrossreff_"+haloname+".csv")                                                        
