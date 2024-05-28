from matplotlib import pyplot as plt
from matplotlib import colors
import math
import pandas as pd


def xyEventDisplay(df, x_scalefactor, y_scalefactor, g_mass, events=None, savefig=False):
    #Plots the xy locations of the CaloHits for each event, with arrows representing the Rhadron momenta. Calohit energies are scaled
    #by their size and color.
    if events is None:
        maxrange = range(max(df['Event']))
    else:
        maxrange = events

    for i in maxrange:
        if events is None:
            j = i+1
        else:
            j = i
        event = df[df['Event'] == j]
        Rhad1_px = x_scalefactor * event['Rhad1_px [GeV]'].iloc[0]
        Rhad1_py = y_scalefactor * event['Rhad1_py [GeV]'].iloc[0]
        Rhad1_ET = (g_mass**2 + (Rhad1_px/x_scalefactor)**2 + (Rhad1_py/y_scalefactor)**2)**0.5
        Rhad2_px = x_scalefactor * event['Rhad2_px [GeV]'].iloc[0]
        Rhad2_py = y_scalefactor * event['Rhad2_py [GeV]'].iloc[0]
        Rhad2_ET = (g_mass**2 + (Rhad2_px/x_scalefactor)**2 + (Rhad2_py/y_scalefactor)**2)**0.5


        plt.clf()
        plt.scatter(event['Calohit X [cm]'], event['Calohit Y [cm]'], s=event['Calohit Energy [GeV]'], 
                    c=event['Calohit Energy [GeV]'], alpha=0.5, cmap='viridis', norm=colors.LogNorm(vmin=0.1, vmax=10000))
        plt.arrow(0,0,Rhad1_px,Rhad1_py, width=1, label=r'Rhadron 1; $E_T={} GeV$'.format(round(Rhad1_ET)), color='orange')
        plt.arrow(0,0,Rhad2_px,Rhad2_py, width=1, label=r'Rhadron 2; $E_T={} GeV$'.format(round(Rhad2_ET)), color='deepskyblue')
        plt.xlabel('X [cm]')
        plt.ylabel('Y [cm]')
        plt.suptitle('Event ' + str(j) + ': CaloHit Locations & Energies for g=1800GeV')
        plt.legend(loc='upper left')
        cbar = plt.colorbar()
        cbar.set_label('Energy [GeV]')
        if savefig:
            plt.savefig('xyEvent' + str(j) + '.png')
        else:
            plt.show()


def rzEventDisplay(df, g_mass, energyCut=0, savefig=False):
    #Plots the rz locations of the CaloHits in the entire dataset
    #reduced_df = df[(df['Calohit Energy [GeV]'] > energyCut) & (df['ECal Type'] == 'EB')]
    plt.clf()
    plt.scatter(df['Calohit Z [cm]'], df['Calohit R [cm]'], alpha=0.5)
    plt.xlabel('Z [cm]')
    plt.ylabel('R [cm]')
    if energyCut == 0:
        plt.suptitle('RZ CaloHit Locations for g=1800GeV')
    else:
        plt.suptitle('RZ CaloHit Locations with Energy > ' + str(energyCut) + ' GeV for g=1800GeV')

    if savefig:
        plt.savefig('rzEventDisplay.png')
    else:
        plt.show()


def hitLocations(df, energyCut=0):
    #Returns the number of hits in the EB, EE, and ES regions for hits with energy > energyCut.
    reduced_df = df[df['Calohit Energy [GeV]'] > energyCut]

    return len(reduced_df[reduced_df['ECal Type']=='EB']), len(reduced_df[reduced_df['ECal Type']=='EE']), len(reduced_df[reduced_df['ECal Type']=='ES'])


def nHitsInEvent(df, energyCut=0):
    #Returns the number of hits in each event with energy > energyCut.
    nHits = []
    for i in range(max(df['Event'])):
        event = df[df['Event'] == i+1]
        nHits.append(len(event[event['Calohit Energy [GeV]'] > energyCut]))
        if len(event[event['Calohit Energy [GeV]'] > energyCut]) == 4:
            print(i+1)
    return nHits


def nHitsAssociatedWithRHadron(df):
    #Returns the number of hits associated with the Rhadrons in the dataset.
    nHits = 0
    nHitsFromRHadron = 0
    for i in range(max(df['Event'])):
        event = df[df['Event'] == i+1]
        Rhad1_px = event['Rhad1_px [GeV]'].iloc[0]
        Rhad1_py = event['Rhad1_py [GeV]'].iloc[0]
        Rhad1_theta = math.atan2(Rhad1_py, Rhad1_px)
        Rhad2_px = event['Rhad2_px [GeV]'].iloc[0]
        Rhad2_py = event['Rhad2_py [GeV]'].iloc[0]
        Rhad2_theta = math.atan2(Rhad2_py, Rhad2_px)
        
        hits = event[(event['Calohit Energy [GeV]'] > 1000)]
        for index, row in hits.iterrows():
            hit_theta = math.atan2(row['Calohit Y [cm]'], row['Calohit X [cm]'])
            if abs(hit_theta - Rhad1_theta) < 0.1 or abs(hit_theta - Rhad2_theta) < 0.1:
                nHitsFromRHadron += 1
            nHits += 1
        
    return nHits, nHitsFromRHadron


def removeMuonHits(df):
    #Removes hits from the muon chamber in the dataset.
    return df[~df['Detector Type'].isin(['MuonRPC', 'MuonDT', 'MuonCSC'])]


def removeECALHits(df):
    #Removes hits from the ECAL in the dataset.
    return df[~df['Detector Type'].isin(['EB', 'EE', 'ES'])]


def analyzeVertices(df):
    spiked_df = df[df['Calo Hit Energy [GeV]'] >= 1000]
    uniqueHits = spiked_df['Calo Hit ID'].unique()
    i=0
    for hit in uniqueHits:
        hit_df = spiked_df[spiked_df['Calo Hit ID'] == hit]
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
            print(hit_df)
        i+=1
        if i == 5:
            break
    print("Number of spiked hits in barrel = {}".format(len(uniqueHits)))


        


file = 'Gluino1800GeV_SimVertices2.csv'
df = pd.read_csv(file)
#x_scalefactor, y_scalefactor = 125 / max(max(df['Rhad1_px [GeV]']), max(df['Rhad2_px [GeV]'])), 125 / max(max(df['Rhad1_py [GeV]']), max(df['Rhad2_py [GeV]']))
#g_mass = 1800

#df = removeMuonHits(df)
#df = removeECALHits(df)
#nHits, nHitsFromRHadron = nHitsAssociatedWithRHadron(df)
#print(nHits, nHitsFromRHadron)
analyzeVertices(df)

#xyEventDisplay(muonless_df, x_scalefactor, y_scalefactor, g_mass, events=None, savefig=False)
#rzEventDisplay(removeMuonHits(df), g_mass, energyCut=0, savefig=False)