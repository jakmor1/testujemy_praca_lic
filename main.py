
import tkinter as tk
import tkinter.filedialog as fd
import pandas as pd
import numpy as np
from IPython.display import display

import matplotlib.pyplot as plt

def load_file(filepath):
    str = filepath.split("/")
    str = str[-1]
    str = str.split(".")
    str = str[1]
    #print(str)


    with open(filepath) as my_file:
        lst = []
        for line in my_file:
            lst.append(line.split())
            #print(line.split())
        #print(lst)
        lst = [[int(x) if x.isdigit() else x for x in subarray] for subarray in lst]  # zamiana string na liczby
        #print(lst)
        df_datafile = pd.DataFrame(data=lst, columns=["Col0", "Col1", "Col2", "Col3", "Col4"])
        #display(df_datafile)
        #print("tutaj")
        #print(df_datafile.loc[1, :].values[0])  # odwołanie do danego wiersza a potem kolumn, a nastepnie do konkretnej wartosci z wiersza

        result = df_datafile[df_datafile["Col1"] == 2]  # wybierz wiersze z danej kolumny, ktore posiadaja dana wartosc
        #display(result)
        #print("pokaz ilosc")
        #print(len(result.index))

        dane = []
        for chromosome in range(1, 6):
            result = df_datafile[df_datafile["Col1"] == chromosome]
            if (len(result.index > 1)):

                starts = result["Col3"].tolist()[:-1]  # bez ostatniego miejsca zakonczenia locus
                ends = result["Col2"].tolist()[1:]  # bez pierwszego rozpoczenia miejsca locus
                num = int(str)

                for i in range(len(starts)):
                    COs = starts[i] + ((ends[i] - starts[i]) // 2)
                    size = ends[i] - starts[i]
                    tmp_lst = [num, chromosome, starts[i], ends[i], COs, size]
                    dane.append(tmp_lst)
        #print(dane)
        df_perfile = pd.DataFrame(data=dane, columns=["Num", "Chr", "Start", "Stop", "Crossover_sites", "Size"])
        #display(df_perfile)
        return df_perfile

def remove_introgression_chr1(df): #remove rows with below condition and reindex
    df = df.drop(df[(df.Chr == 1) & (df.Crossover_sites > 9018718) & (df.Crossover_sites < 10012987)].index)
    df.reset_index(inplace=True, drop=True)
    return df

def occurrences(df): #count crossing-over occurences for every unique lib(column "Num")
    lst = df["Num"].tolist()
    lst = list(set(lst))

    occurrences = []
    for unique_num in lst:
        count = df['Num'].value_counts()[unique_num]
        occurrences.append(count)
    return occurrences







root = tk.Tk()
root.withdraw()


control_size = 0
mutant_size = 0


file_paths_control = fd.askopenfilenames(parent=root, title="Kontrola")

control_COs = pd.DataFrame(columns=["Num","Chr","Start","Stop","Crossover_sites","Size"])
for filepath in file_paths_control:
    ret = load_file(filepath)
    frames = [control_COs, ret]
    control_COs = pd.concat(frames,ignore_index=True)
    control_size +=1
    #print("pokaz mi")



file_paths_mutant = fd.askopenfilenames(parent=root, title="Mutant")

mutant_COs = pd.DataFrame(columns=["Num","Chr","Start","Stop","Crossover_sites","Size"])
for filepath in file_paths_mutant:
    ret = load_file(filepath)
    frames = [mutant_COs, ret]
    mutant_COs = pd.concat(frames,ignore_index=True)
    mutant_size +=1
    #print("pokaz mi")

control_COs.to_csv("contr_all.csv")  #nwm czy potrzebne sa te zapisy
mutant_COs.to_csv("mut_all.csv")

display(control_COs)
display(mutant_COs)

control_COs = remove_introgression_chr1(control_COs)
mutant_COs = remove_introgression_chr1(mutant_COs)

print("po usunieciu introgression")
display(control_COs)
display(mutant_COs)

print(len(file_paths_control))
print(len(file_paths_mutant))


occurrences_control = occurrences(control_COs)
occurrences_mutant = occurrences(mutant_COs)

print(occurrences_control)
print(occurrences_mutant)

print(len(occurrences_control))
print(len(occurrences_mutant))


from scipy import stats
# przeprowadzenie t-testu
# t_stat, p_val = stats.ttest_ind(occurrences_mutant, occurrences_control)
#
# # wyświetlenie wyników
# print("T-stat: ", t_stat)
# print("P-value: ", p_val)
q1, q3 = np.percentile(occurrences_control, [25, 75])
iqr = q3 - q1
lower_bound = q1 - 1.5*iqr
upper_bound = q3 + 1.5*iqr
occurrences_control = [x for x in occurrences_control if lower_bound <= x <= upper_bound]

q1, q3 = np.percentile(occurrences_mutant, [25, 75])
iqr = q3 - q1
lower_bound = q1 - 1.5*iqr
upper_bound = q3 + 1.5*iqr
occurrences_mutant = [x for x in occurrences_mutant if lower_bound <= x <= upper_bound]


counts_control = dict()
for i in occurrences_control:
  counts_control[i] = counts_control.get(i, 0) + 1

counts_mutant = dict()
for i in occurrences_mutant:
  counts_mutant[i] = counts_mutant.get(i, 0) + 1

print(counts_control)
print(counts_mutant)



#jeszcze moze dodac warunek if w zaleznosci co wieksze bedzie by mniejsze na wierzchu bylo

plt.bar(counts_control.keys(), counts_control.values(),  color='g',align='center', label = "control")
plt.bar(counts_mutant.keys(), counts_mutant.values(),color='r',align='center', label = "mutant")
plt.xlabel('Crossovers',fontsize=15)
plt.ylabel('No of plants',fontsize=15)
plt.legend(loc="upper right")
plt.title(label='Crossover frequency',fontsize=25)
plt.savefig("Crossover frequency")
plt.close()


centromeres = [15086045,3607929,13587786,3956021,11725024]
chrs = ["Chr1","Chr2","Chr3","Chr4","Chr5"]
chr_ends = [30427671,19698289,23459830,18585056,26975502]

#control_size,mutant_size


#df = df.drop(df[(df.Chr == 1) & (df.Crossover_sites > 9018718) & (df.Crossover_sites < 10012987)].index)
from scipy import signal
ma = [0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125]
for chromosome in range(1,6): #dla kazdego chromosomu
    tmp_control = control_COs[control_COs.Chr == chromosome] #wez wiersze tylko z danym nr chromosomu
    tmp_mutant = mutant_COs[mutant_COs.Chr == chromosome]
    window = list(range(1, chr_ends[chromosome-1], 300000)) #stworz okna z krokiem 300000 oraz dodaj koniec chromosomu do okna
    window.append(chr_ends[chromosome-1])

    last_el_in_window = (chr_ends[chromosome-1] - window[len(window)-2]) /300000

    print("pprawda")
    print(last_el_in_window)

    window_control = []
    window_mutant = []
    for i in range(len(window)): #przejedz po kazdej wartosci z okna tj. 1,300001, 600001,... aż do końca i podczas każdego przejechania policz ilosc COs w danym oknie
        #print(i)
        #print(window[i])
        COs_in_window = len(tmp_control[(tmp_control.Crossover_sites >= window[i-1]) & (tmp_control.Crossover_sites <= window[i])])
        window_control.append(COs_in_window)
        COs_in_window = len(tmp_mutant[(tmp_mutant.Crossover_sites >= window[i-1]) & (tmp_mutant.Crossover_sites <= window[i])])
        window_mutant.append(COs_in_window)

    print("mmmooj")
    print(window_control[len(window_control)-1])
    window_control[len(window_control)-1] =  window_control[len(window_control)-1]/last_el_in_window
    print(window_control[len(window_control) - 1])
    window_mutant[len(window_mutant) - 1] = window_mutant[len(window_mutant) - 1] / last_el_in_window
    window_control = [x / mutant_size for x in window_control]
    window_mutant = [x / mutant_size for x in window_mutant]

    # window_control.insert(0,0)
    # window_control.insert(0, 0)
    # window_control.insert(0, 0)
    # window_control.append(0)
    # window_control.append(0)
    # window_control.append(0)
    #
    # window_mutant.insert(0,0)
    # window_mutant.insert(0, 0)
    # window_mutant.insert(0, 0)
    # window_mutant.append(0)
    # window_mutant.append(0)
    # window_mutant.append(0)

    filt_control = signal.convolve(window_control, ma, mode='same')
    filt_mutant = signal.convolve(window_mutant, ma, mode='same')

    #filt_control =window_control
    #filt_mutant = window_mutant
    # filt_control = filt_control[3:-3]
    # filt_mutant = filt_mutant[3:-3]


    print(type(filt_control))
    filt_control = np.multiply(filt_control,(186 / 1859))
    print(type(filt_control))

    print(filt_control)
    print(filt_mutant)
    plt.plot(window, filt_control, color='r', label="control")
    plt.plot(window, filt_mutant, color='b', label="mutant")
    plt.axvline(centromeres[chromosome-1], linestyle='--', color='g', label="centromere")

    plt.title("COs distribution for Chr" + str(chromosome))
    plt.legend(loc="upper right")
    plt.savefig("COs distribution for Chr" + str(chromosome))
    plt.close()

    print(filt_control)
    print(filt_mutant)



#
# display(tmp_control)
# display(tmp_mutant)
# print(window)
# print(len(window))
#
# print(window_control)
# print(len(window_control))
#
# print(window_mutant)
# print(len(window_mutant))
#
# print(filt_control)
# print(filt_mutant)
#
# print(len(filt_control))
# print(len(filt_mutant))

# plt.plot(window, filt_control,  color='r',label ="control")
# plt.plot(window, filt_mutant,color='b',label="mutant")
# plt.axvline(centromeres[4],linestyle='--', color='g', label = "centromere")
#
# plt.title("COs distribution for Chr" )
# plt.legend(loc="upper right")
# plt.savefig("COs distribution for Chr")
# plt.show()




#
#
# df = pd.DataFrame(data = [[1,4,105,178,145,40]],columns=["Num", "Chr", "Start","Stop","Crossover_sites","Width"])
# display(df)
#
# print()
# print("tutaj")
# for file_path in file_paths_control:
#     print(file_path)
#
# for file_path in file_paths_mutant:
#     print(file_path)
#
#
#
#
#
#
#
# for i in range(1,6):
#     print(i)
#









        # df_per_chr = pd.DataFrame(data=lst, columns=["Num","Chr","Start","Stop","Crossover_sites","Size"])

        # result = result.append([12,12,12,12,"XD"])  ignore index i moze jeszcze reset index
        # display(result)
        # starts = result["Col3"].tolist()
        # ends = result["Col2"].tolist()
        # print(starts)
        # print(ends)


# #str = "C:/Users/Kuba/Downloads/licencjat/test/rad17.108.mask.smooth.co.txt"
# str = "C:/Users/Kuba/Downloads/licencjat/test/beth.3484102.mask.smooth.co.txt"
# str = str.split("/")
# str = str[-1]
# str = str.split(".")
# str = str[1]
# print(str)
#
#
#
#
#
#
# #with open("C:/Users/Kuba/Downloads/licencjat/test/rad17.108.mask.smooth.co.txt") as my_file:
# with open("C:/Users/Kuba/Downloads/licencjat/test/beth.3484102.mask.smooth.co.txt") as my_file:
#     lst=[]
#     for line in my_file:
#         lst.append(line.split())
#         print(line.split())
#     print(lst)
#     lst = [[int(x) if x.isdigit() else x for x in subarray] for subarray in lst]  #zamiana string na liczby
#     print(lst)
#     df_datafile = pd.DataFrame(data=lst,columns=["Col0", "Col1", "Col2", "Col3", "Col4"])
#     display(df_datafile)
#     print("tutaj")
#     print(df_datafile.loc[1,:].values[0])   #odwołanie do danego wiersza a potem kolumn, a nastepnie do konkretnej wartosci z wiersza
#
#     result = df_datafile[df_datafile["Col1"]==2] #wybierz wiersze z danej kolumny, ktore posiadaja dana wartosc
#     display(result)
#     print("pokaz ilosc")
#     print(len(result.index))
#
#     dane = []
#     for chromosome in range(1,6):
#         result = df_datafile[df_datafile["Col1"] == chromosome]
#         if(len(result.index>1)):
#
#
#             starts = result["Col3"].tolist()[:-1] #bez ostatniego miejsca zakonczenia genu
#             ends = result["Col2"].tolist()[1:] #bez pierwszego rozpoczenia miejsca genu
#             num = int(str)
#
#
#             for i in range(len(starts)):
#                 COs =starts[i] + ((ends[i] - starts[i])//2)
#                 size = ends[i]-starts[i]
#                 tmp_lst = [num,chromosome,starts[i],ends[i],COs,size]
#                 dane.append(tmp_lst)
#     print(dane)
#     df_perfile = pd.DataFrame(data=dane,columns=["Num","Chr","Start","Stop","Crossover_sites","Size"])
#     display(df_perfile)
#
#
#
#            # df_per_chr = pd.DataFrame(data=lst, columns=["Num","Chr","Start","Stop","Crossover_sites","Size"])
#
#
#
#
#
#     #result = result.append([12,12,12,12,"XD"])  ignore index i moze jeszcze reset index
#     display(result)
#     starts = result["Col3"].tolist()
#     ends = result["Col2"].tolist()
#     print(starts)
#     print(ends)





















# import numpy as np
# import matplotlib.pyplot as plt
#
#
# fig, (ax1, ax2) = plt.subplots(2, 1)
# # make a little extra space between the subplots
# fig.subplots_adjust(hspace=0.5)
#
# dt = 0.01
# t = np.arange(0, 30, dt)
#
# # Fixing random state for reproducibility
# np.random.seed(19680801)
#
#
# nse1 = np.random.randn(len(t))                 # white noise 1
# nse2 = np.random.randn(len(t))                 # white noise 2
# r = np.exp(-t / 0.05)
#
# cnse1 = np.convolve(nse1, r, mode='same') * dt   # colored noise 1
# cnse2 = np.convolve(nse2, r, mode='same') * dt   # colored noise 2
#
# # two signals with a coherent part and a random part
# s1 = 0.01 * np.sin(2 * np.pi * 10 * t) + cnse1
# s2 = 0.01 * np.sin(2 * np.pi * 10 * t) + cnse2
#
# ax1.plot(t, s1, t, s2)
# ax1.set_xlim(0, 5)
# ax1.set_xlabel('Time')
# ax1.set_ylabel('s1 and s2')
# ax1.grid(True)
#
# cxy, f = ax2.csd(s1, s2, 256, 1. / dt)
# ax2.set_ylabel('CSD (dB)')
# plt.show()


