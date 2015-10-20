#!/env python

import os.path, re, string
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

####################################################################################################
### Variable definitions
####################################################################################################

clust_header = ("Volume", "CM_RL", "CM_AP", "CM_IS", "minRL",
                "maxRL", "minAP", "maxAP", "minIS", "maxIS", "Mean", "SEM", "Max Int",
                "MI_RL", "MI_AP", "MI_IS")

task="restingstate"
usedFwhm="4.2"

root_dir="/data"
data_dir=os.path.join(root_dir, "sanDiego", "rsfcGraphAnalysis", "data")
admin_data_dir=os.path.join(data_dir, "admin")
config_data_dir=os.path.join(data_dir, "config")
seeds_data_dir=os.path.join(data_dir, "seeds")

groups="mddOnly"

####################################################################################################
### Function definitions
####################################################################################################


def substituteShortLabels(inLevel):
    returnSubstitutedLabel = re.sub("[0-9]+ ", "", re.sub("Inf", "Inferior", re.sub("Sup", "Superior", re.sub("Gy", "Gyrus",
        re.sub("^R", "Right", re.sub("^L", "Left", inLevel))))))
    
    return (returnSubstitutedLabel)



def makePublicationTable(inClusterWhereAmI, inClusters, inRoistats, inRoistats_averageStatValue=None, inRoistats_averageCoefficientValue=None,
                         inStatColumnName="Default Stat Name", inCoefficientColumnName="Default Coefficient Name", inCom=True):

    def extractHemisphereLabels(xx):
        return(re.sub("[^RL]", "", xx[0]))

    def stripHemisphereLabels(xx):
        return(re.sub("[RL] ", "", xx))

    hemisphere=map(extractHemisphereLabels, inClusterWhereAmI)
    ## print hemisphere

    if inCom:
        locations=pd.concat([pd.DataFrame(map(stripHemisphereLabels, inClusterWhereAmI)),
                             pd.DataFrame(hemisphere),
                             np.round(pd.DataFrame(inClusters[["Volume", "CM_RL", "CM_AP", "CM_IS"]]), 0)],
                            axis=1)
        locations.columns=["Structure", "Hemisphere", "Volume", "CM_RL", "CM_AP", "CM_IS"]
    else:
        locations=pd.concat([pd.DataFrame(map(stripHemisphereLabels, inClusterWhereAmI)),
                             pd.DataFrame(hemisphere),
                             np.round(pd.DataFrame(inClusters[["Volume", "MI_RL", "MI_AP", "MI_IS"]]), 0)],
                            axis=1)
        locations.columns=["Structure", "Hemisphere", "Volume", "MI_RL", "MI_AP", "MI_IS"]

    r=re.compile("^Mean_")
    meanCols=filter(r.match, inRoistats.columns)

    if len(meanCols) == 0:
        raise Error("There are no columns in the inRoistats data frame that begin with Mean_")
    
    inRoistats_grouped = inRoistats.groupby("Group")
    means=inRoistats_grouped[meanCols].mean().transpose().reset_index(drop=True)

    locations = pd.concat([locations, means], axis=1)

    if inRoistats_averageStatValue is not None:
        # drop the sub-brik (0'th) column from the stat data frame
        inRoistats_averageStatValue.drop(inRoistats_averageStatValue.columns[[0]], axis=1, inplace=True)        
        locations=pd.concat([locations, inRoistats_averageStatValue.transpose().reset_index(drop=True)], axis=1)
        locations.rename(columns={0:inStatColumnName}, inplace=True)

    if inRoistats_averageCoefficientValue is not None:
        # drop the sub-brik (0'th) column from the coefficient data frame        
        inRoistats_averageCoefficientValue.drop(inRoistats_averageCoefficientValue.columns[[0]], axis=1, inplace=True)            
        locations=pd.concat([locations, inRoistats_averageCoefficientValue.transpose().reset_index(drop=True)], axis=1)
        locations.rename(columns={0:inCoefficientColumnName}, inplace=True)        

    return locations
                                
def savePublicationTable(inPublicationTable, inPublicationTableFilename, append=True):
    if append:
        print "*** Appending publication table to", inPublicationTableFilename        
        inPublicationTable.to_csv(inPublicationTableFilename, sep=",", index=False, mode="a")
    else:
        print "*** Writing publication table to", inPublicationTableFilename                
        inPublicationTable.to_csv(inPublicationTableFilename, sep=",",  index=False)

def readStatsTable(inFilename):

    print ("*** Reading %s" % inFilename)
    statsTable=pd.read_table(inFilename, comment="#", sep=r"\s*")
    ## dump the first column as it's only the file name
    ## statsTable=statsTable[, -1]
    statsTable.drop(statsTable.columns[[0]], axis=1, inplace=True)

    return(statsTable)


def readClustersTable(inFilename):
    print("*** Reading %s" % (inFilename))
    clusters=pd.read_table(inFilename, comment="#", sep="\s*", names=clust_header)

    return (clusters)


def readClusterLocationsTable(inFilename):
    print "*** Reading cluster locations from %s" % inFilename
    ## the gsub here chews up multiple consequtive spaces and replaces them with a single space
    ff=open(inFilename, "r")
    ll=string.strip(ff.readline())
    ff.close()
    # print "ll: ", ll
    clusterWhereAmI=re.split(",\s*", ll)
    # print (clusterWhereAmI)
    clusterWhereAmI=map(lambda x: re.sub(" +", " ", x), clusterWhereAmI)
    #   clusterWhereAmI=gsub(" +", " ", scan(file=inFilename, what='character', sep=',', quiet=TRUE))
    
    return (clusterWhereAmI)


def readSubjectOrderTable(inFilename):
    print("*** Reading %s" % inFilename)
    subjectOrder=pd.read_csv(inFilename, comment="#")

    return(subjectOrder)


def fixSubjectOrderTable(inSubjectOrderTable):
    inSubjectOrderTable['subject']=map(lambda x: re.sub("_A", "", x), inSubjectOrderTable['subject'])
    
    ## inSubjectOrderTable[inSubjectOrderTabl['$subject']=="300", "subject"] = "169/300"
    ## inSubjectOrderTable$subject=as.factor(inSubjectOrderTable$subject)

    return(inSubjectOrderTable)


## Reads the seed file and does the (crude) equivalent of BAS variable
## substitution
def readSeedsFile(inSeedsFile):
    if os.path.isfile(inSeedsFile):
        print("*** Reading seed from %s\n" % (inSeedsFile))
        ff=open (inSeedsFile, "r")
        table=ff.readlines()
        table=map(string.strip, table)
        print type(table)
        ff.close()
        
        def _subData(inLine):
            return re.sub("$DATA", data_dir, inLine)
        
        table=map(_subData, table)
    else:
        print ("*** readSeedsFile No such file: %s" % (inSeedsFile))
        table=None
               
    return (table)


## extracts the seed name from a file path name pointing to a NIfTI
## file containing the seed
def getSeedName(inSeedPath):
    name=os.path.basename(inSeedPath)
    if (re.search("\\.nii", name)):
        name=re.sub("\\.nii.*", "", name)
    elif (re.search("\\+tlrc", name)):
        name=re.sub("\\+tlrc.*", "", name)
    else:
        pass
    
    return (string.strip(name))

def readCsvFile(inFilename, inSubjectColumnName="ID"):

    print "*** Reading", inFilename
    rCsv=pd.read_csv(inFilename, comment="#", na_values = ["NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT", ".", ""])
    print "*** Read data for %s rows" % (len(rCsv[inSubjectColumnName]))

    return(rCsv)


def read_change_score_file(in_variable):

    change_score_filename=os.path.join(admin_data_dir, "new.mdd.%s.change.score.csv" % in_variable)
    print "*** Change score file is %s" % change_score_filename
    if (os.path.isfile(change_score_filename)):
        change_scores=readCsvFile(change_score_filename)
    else:
        raise IOError("Couldn't find the change score file to go with a variable named " % in_variable)

    return(change_scores)


def getClusterCount(in_columns_names):
    r=re.compile("^Mean_")
    meanCols=filter(r.match, in_columns_names)
   
    return(len(meanCols))

def generateGraphs(group_data_dir, group_results_dir, rvVariable, rvName, publicationTableFilename, seed_list):

    if groups == "mddAndCtrl" and re.search("diff", rvVariable):
        print("*** You are trying to regress/graph the control and MDD subjects aagainst a time A to C difference. You cannot do that!")
        return()

    print seed_list
    
    for seed in seed_list:
    #for seed in [seed_list[0]]:        
        seedName=getSeedName(seed)

        print "SeedName: %s" % seedName

        try:
            ff=open(publicationTableFilename, "w+")
            ff.write("Seed=%s,Variable=%s (%s),\n" % (seedName, re.sub("\n", " ", rvName), rvVariable))
        except Error as (errno, strerror):
            print "*** Couldn't write to publication file: {0}: {1}".format(errno, strerror)
            raise 
        else:
            ff.close()
            
        infix="regression.fwhm%s.%s.%s.%s.and.%s" % (usedFwhm, task, groups, seedName, rvVariable)

        roistats_filename=os.path.join(group_results_dir, "roiStats.%s.txt" % (infix) )
        roistats_averageTvalue_filename=os.path.join(group_results_dir, "roiStats.%s.averageTValue.txt" % (infix) )
        roistats_averageCoefficientValue_filename=os.path.join(group_results_dir, "roiStats.%s.averageCoefficientValue.txt" % (infix))

        if os.path.isfile(roistats_filename):
            roistats=readStatsTable(roistats_filename)
            roistats_averageTvalue=readStatsTable(roistats_averageTvalue_filename)
            roistats_averageCoefficientValue=readStatsTable(roistats_averageCoefficientValue_filename)

            clusterCount=getClusterCount(roistats.columns)
            if clusterCount > 0:
                print "*** %d ** clusters found in %s" % (clusterCount, roistats_filename)

                clustersFilename=os.path.join(group_results_dir, "clust.%s.txt" % infix)
                clusters=readClustersTable(clustersFilename)

                ## print(clusters)
                
                ## this table contains the locations, as text, of the clusters and is the output of a perl script
                clusterLocationsFilename=os.path.join(group_results_dir, "clusterLocations.%s.csv" % infix)
                clusterWhereAmI=readClusterLocationsTable(clusterLocationsFilename)
                # print clusterWhereAmI
                
                ## this file stores the order of the subjects in each of the following BRIK files
                subjectOrderFilename=os.path.join(group_data_dir, "subjectOrder.%s.%s.csv" % (groups, seedName))
                subjectOrder=fixSubjectOrderTable(readSubjectOrderTable(subjectOrderFilename))
                subjectOrder['subject']=subjectOrder['subject'].astype(int)

                print "*** Read subject order data for %s subjects" %  len((subjectOrder.index))

                change_scores=read_change_score_file(rvVariable)
                change_scores['ID']=change_scores['ID'].astype(int)

                mgd=pd.concat([subjectOrder, roistats], axis=1)

                mgd=pd.merge(mgd, change_scores, left_on="subject", right_on="ID", sort=False)
                mgd.drop('ID', axis=1, inplace=True)
                mgd.rename(columns= {"Grp": "Group"}, inplace=True)

                publicationTable=makePublicationTable(clusterWhereAmI, clusters, mgd, roistats_averageTvalue, roistats_averageCoefficientValue,
                                                      inStatColumnName="Average t value", inCoefficientColumnName="Average Coefficient Value", inCom=True)
                savePublicationTable(publicationTable, publicationTableFilename, True)
                
                ## print(publicationTable)
                
                melted_mgd=pd.melt(mgd,
                                   id_vars=["subject", "Group", rvVariable],
                                   value_vars=["Mean_" + str(ii) for ii in range(1, clusterCount+1)],
                                   var_name="cluster")

                uniqueClusterNames= ", ".join("%02d %s" % t for t in zip (range(1, len(clusterWhereAmI)+1), clusterWhereAmI)).split(", ")
                clusterToTalairachNames=dict(zip(["Mean_" + str(ii) for ii in range(1, len(uniqueClusterNames)+1)], uniqueClusterNames))
                # print clusterToTalairachNames
                # quite equivalent to paste in R
                # 
                
                melted_mgd['cluster'] =melted_mgd['cluster'].astype('category')
                #melted_mgd['cluster'] = pd.Category[
                
                ## print melted_mgd
                ## print melted_mgd.info()

                graphRegressions(melted_mgd, group_results_dir, rvVariable, rvName, seedName, clusterToTalairachNames)
        else:
            print "*** No such file: %s" % roistats_filename    

def graphRegressions(melted_mgd, group_results_dir, rvVariable, rvName, seedName, clusterToTalairachNames):

    ## print(head(melted.mgd))
    imageDirectory=os.path.join(group_results_dir, seedName)
    if not os.path.isdir(imageDirectory):
        os.makedir(imageDirectory, 0770)

    ## print "Categories =>", [melted_mgd['cluster'].cat.categories[0]]

    for level in melted_mgd['cluster'].cat.categories:
    #for level in [melted_mgd['cluster'].cat.categories[0]]:        

        ss=melted_mgd[melted_mgd['cluster'] == level]
        # print ss
        
        talairachName=clusterToTalairachNames[level]
        # print "%s => %s" % (level, talairachName)
        imageFilename=os.path.join(imageDirectory, "%s.fwhm%s.%s.py.pdf" % (re.sub(" +", ".", talairachName),  usedFwhm, rvVariable))
        print "*** Creating %s" % imageFilename

        roistats_summary=ss[rvVariable].describe()
        ## print roistats_summary
        if re.search("\\.reversed$", group_results_dir):
            print "*** Using reversed axis variables and variable names"
            x_axis="value"
            y_axis=rvVariable
            x_axis_label="RSFC (Z-score)"
            y_axis_label=rvName
        else:
            x_axis=rvVariable
            y_axis="value"
            x_axis_label=rvName
            y_axis_label="RSFC (Z-score)"

        figure = plt.figure()
        ax = figure.add_subplot(111)
        ax.set_title(substituteShortLabels(talairachName))
        ax.set_xlabel(x_axis_label)
        ax.set_ylabel(y_axis_label)
                     
        plt.scatter(ss[x_axis], ss[y_axis], c="red")
        plt.hold(True)
        # fit an order 1 polynomial (linear) model to the the x data
        fit = np.polyfit(ss[x_axis], ss[y_axis], 1)
        # create a function (fit_fn) that returns the fit values for x
        # based on the fit model in fit
        fit_fn = np.poly1d(fit)
        plt.plot(ss[x_axis], fit_fn(ss[x_axis]), "black")
        plt.show()
        #plt.savefig(imageFilename)
        
    
def main():

    regressionVariables=[
        # ("CDRS.t.score.diff",             "Children's Depression Rating Scale\n(Baseline to 3 Months Change)"),
        # ("MASC.tscore.diff",              "Multidimensional Anxiety Scale for Children\n(Standardized; Baseline to 3 Months Change)"),
        # ("CGAS.diff",                     "Children's Global Assessment Scale\n(Baseline to 3 Months Change)"),
        # ("RADS.Total.tscore.diff",        "Reynolds Adolescent Depression Scale Total\n(Standardized; Baseline to 3 Months Change)"),
        # ("CDRS.t.score.scaled.diff",      "Children's Depression\nRating Scale\n(Baseline to 3 Months Change)")#,
        
        ("CDRS.t.score.scaled.diff",         r"$\Delta$ CDRS-R")#,
        
        # ("CGAS.scaled.diff",              "Children's Global Assessment Scale\n(Baseline to 3 Months Change)")
        # ("MASC.tscore.scaled.diff",       "Multidimensional Anxiety Scale for Children\n(Standardized; Baseline to 3 Months Change)"),
        # ("RADS.Total.tscore.scaled.diff", "Reynolds Adolescent Depression Scale Total\n(Standardized; Baseline to 3 Months Change)")
        # ("BDI.II.Total.diff",             "Beck Depression Inventory II (A to C Change)"),
        # ("CDI.Total.diff",                "Children's Depression Inventory (A to C Change)")
        ]

    #seedFile=os.path.join(config_data_dir, "juelich_whole_amygdala_seeds.txt")
    seedFiles=[os.path.join(config_data_dir, "juelich_amygdala_seeds_weights.txt"), ]

    print seedFiles
    for seedFile in seedFiles:
        seed_list=readSeedsFile(seedFile)

        for rvVariable, rvName in regressionVariables:
        
            group_data_dir=os.path.join(data_dir, "Group.data.%s.withAandC" % (rvVariable))
            group_results_dir=os.path.join(data_dir, "Group.results.%s.withAandC.reversed" % (rvVariable))

            ## publicationTableFilename=os.path.join(group_results_dir, "publicationTable.%s.%s.csv.py" % (usedFwhm, groups))        
            publicationTableFilename="publicationTable.%s.%s.csv.py" % (usedFwhm, groups)

            generateGraphs(group_data_dir, group_results_dir, rvVariable, rvName, publicationTableFilename, seed_list)

