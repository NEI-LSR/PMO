# Written with the help of ChatGPT

import pandas as pd
import pickle

# Load the dictionary from the .pkl file
with open('measurements/data_init_1.pkl', 'rb') as file:
    data_dict = pickle.load(file)
df_power = pd.DataFrame()

# Populating the DataFrame
for key, values in data_dict.items():
    df_power[f'Key_{key}'] = values['power']

# Save the DataFrame to a new CSV file
csv_filepath_power = 'measurements/data_init_1.csv'

df_power.to_csv(csv_filepath_power, index=False)
