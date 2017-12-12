import CournotMerge as cm

# Define the markets (use the same order for elasticities, prices, and firm quantities below)
markets = ["San Diego", "Seattle", "DC"]

# Define the market elasticities
elasticity = [0.1, 0.6, 1.1]

# Define the market prices
price = [0.66, 1.20, 1.90]

# Define the Firms as a dictionary starting with index 1 (Name, quantities, markets)
firms={}
firms[1] = cm.firm("Trader Joe's", [45, 25, 15], markets)
firms[2] = cm.firm("Albertson's", [40, 18, 6], markets)
firms[3] = cm.firm("Whole Foods", [55, 45, 19], markets)
firms[4] = cm.firm("Sprouts", [11, 7, 4], markets)

# Identify the indices for the two parties that are merging
merge_index = [1,3]

# Identify any MC efficiencies for the merging parties (efficiency% reduction in MC)
efficiency = 10

# Call the merger program
CM = cm.cournotMerge(F=firms, markets=markets, e=elasticity, p=price, MF=merge_index, E=efficiency)

# Describe the merger simulation
print(CM.describe())

# Summarize the merger sim
CM.summarize()
