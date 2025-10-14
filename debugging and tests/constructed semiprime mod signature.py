import pandas as pd
import json

# Load your data (replace 'your_file.xlsx' with the actual file path)
df = pd.read_excel("mod 2310.xlsx", index_col=0)  # Assuming the first column is the index

# Step 1: Clean the rows where row header is divisible by 5
df = df[df.index % 5 != 0]

# Step 2: Clean the columns where column header is divisible by 5
df = df.loc[:, df.columns % 5 != 0]

# Now 'df' should be cleaned
print(df.head())

valid_headers = [header for header in df.columns if header % 5 != 0]
valid_indexes = [index for index in df.index if index % 5 != 0]

# Create the dictionary dynamically with valid keys
result_dict = {header: set() for header in valid_headers}  # Use a set to store unique values

# Step 2: Loop through the DataFrame and populate the dictionary
for row_index, row in df.iterrows():  # row_index is the row header
    for col_index, value in row.items():  # col_index is the column header
        if value in result_dict:
            # Compute (column_header + row_index) % 2310 and add to the set
            result_dict[value].add((col_index + row_index) % 2310)

# Convert sets to sorted lists for each key in the dictionary
for key in result_dict:
    result_dict[key] = sorted(result_dict[key])

# Step 3: Print the result_dict to see the output (printing the first 5 items)
for key, value in list(result_dict.items())[:5]:
    print(f"{key}: {value}")
for key, value in list(result_dict.items())[-5:]:
    print(f"{key}: {value}")

with open('2310 signature.json', 'w') as json_file:
    json.dump(result_dict, json_file, indent=4)
with open('2310 signature.txt', 'w') as txt_file:
    txt_file.write(str(result_dict))
