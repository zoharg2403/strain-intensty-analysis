import pandas as pd
import os
import re

# import raw data .txt files
raw_data_folder = r'C:\Users\zoharga\Desktop\Zohar\plot\Raw data'
raw_data_files = os.listdir(raw_data_folder)
# import strain data file
ORF_Gene_data = pd.read_excel("List of genes 210 x Tef Cherry oex Javier.xlsx",
                              usecols=["ORF", "Gene", "384 Plate_Row_Col"])
ORF_Gene_data['Plate'] = ORF_Gene_data["384 Plate_Row_Col"].apply(lambda s: list(re.split('(\d+)', s))[1])

# create list of dataframes for each data file
df_list = []
for f in raw_data_files:
    # get cur_plate data:
    cur_plate_data = pd.read_csv(os.path.join(raw_data_folder, f), sep='\t')
    cur_plate_strain_data = ORF_Gene_data[ORF_Gene_data['Plate'] == str(int(''.join(filter(str.isdigit, f))))]
    cur_plate_data['Plate'] = [str(int(''.join(filter(str.isdigit, f))))] * cur_plate_data.shape[0]

    # remove objects out of gate ('R01')
    cur_plate_data = cur_plate_data[cur_plate_data['R01'] == 1]

    # add 'Gene' and 'ORF' columns to cur_plate_data from ORF_Gene_data
    wells_col = cur_plate_data['Well '].replace(list(range(1, 385)), cur_plate_strain_data["384 Plate_Row_Col"])
    cur_plate_data['Gene'] = wells_col.replace(list(cur_plate_strain_data["384 Plate_Row_Col"]),
                                               list(cur_plate_strain_data['Gene']))
    cur_plate_data['ORF'] = wells_col.replace(list(cur_plate_strain_data["384 Plate_Row_Col"]),
                                              list(cur_plate_strain_data['ORF']))

    # add cur_plate_data to df_list
    df_list.append(cur_plate_data)

# merge a list of dataframes df_list to create one dataframe 'df' (and save to .csv file)
df = pd.concat(df_list, ignore_index=True)
# df.to_csv('full data - plates 1-16.csv', index=False)

# create new dataframe 'strain_mean_intensity' that contains mean intensity 488nm for each strain in 'df'
data_per_strain = pd.DataFrame({
    "ORF": ORF_Gene_data['ORF'],
    "Gene": ORF_Gene_data['Gene'],
    "Well ID": ORF_Gene_data['384 Plate_Row_Col'],
    "Well cell count": ["None"] * len(ORF_Gene_data['ORF']),
    "Mean Intensity 488nm": ["None"] * len(ORF_Gene_data['ORF']),
    "STD Intensity 488nm": ["None"] * len(ORF_Gene_data['ORF']),
    "Median Intensity 488nm": ["None"] * len(ORF_Gene_data['ORF'])})

counter = len(ORF_Gene_data['Gene'])

for G in ORF_Gene_data['Gene']:
    cur_strain_data = df[df['Gene'] == G]
    data_per_strain.loc[data_per_strain['Gene'] == G, "Well cell count"] = cur_strain_data.shape[0]
    data_per_strain.loc[data_per_strain['Gene'] == G, 'Mean Intensity 488nm'] = cur_strain_data["Mean Intensity 488nm"].mean()
    data_per_strain.loc[data_per_strain['Gene'] == G, 'STD Intensity 488nm'] = cur_strain_data["Mean Intensity 488nm"].std()
    data_per_strain.loc[data_per_strain['Gene'] == G, 'Median Intensity 488nm'] = cur_strain_data["Mean Intensity 488nm"].median()
    counter -= 1
    print(counter)

# sort strain_mean_intensity by ascending mean intensity order
data_per_strain.dropna(inplace=True)
data_per_strain = data_per_strain.sort_values(by='Mean Intensity 488nm')
# save strain_mean_intensity to .csv file
data_per_strain.to_csv('Data per Strain.csv', index=False)


