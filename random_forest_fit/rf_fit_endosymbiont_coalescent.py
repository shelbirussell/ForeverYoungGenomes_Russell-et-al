from sklearn.ensemble import RandomForestRegressor
import sklearn.model_selection as xval
import forestci as fci
import numpy as np
import matplotlib 
matplotlib.use('agg') 
import matplotlib.pyplot as plt
import math

# read data and split into variables to predict and summary statistics
data = np.loadtxt("20_stats.txt", delimiter="\t")
#load our species 
stats = np.loadtxt("BathyNewEngland.data", delimiter="\t")
stats = stats[:,[0,2,3,4,5,6,7,8,9,10,11]]

#Data is:
#1. rho
#2. length of rec
#3. theta

#Summaries are:
#4. pi
#5. s
#6. prop 4d < 1kb
#7. 1-3kb
#8. 3-5kb
#9. 5-15kb
#10. 15-30kb
#11. 30-50kb
#12. 50-75kb
#13. 75-100kb
#14. 100-200kb
#15. gt200kb


## log transfor some variables
#data[:,[0]] = np.log10(data[:,[0]])
#data[:,[1]] = np.log10(data[:,[1]])
#data[:,[2]] = np.log10(data[:,[2]])

### get the prediction and data summaries split
data_predict = np.log10(data[:,[2]])
#data_predict = np.log10(data[:,[0]]*data[:,[1]]) 
#data_summary = data[:,[3,4,5,6,7,8,9,10,11,12,13,14]]
data_summary = data[:,[3,5,6,7,8,9,10,11,12,13,14]] 

# split data into training and test set, for now, we'll just keep 10% for testing
data_summary_train, data_summary_test, data_predict_train, data_predict_test = xval.train_test_split(data_summary, data_predict, test_size=0.01, random_state=42)

# create RandomForestRegressor
forest = RandomForestRegressor(n_estimators=1000, random_state=42, oob_score=True, max_depth=8, max_features=10 )
forest.fit(data_summary_train, data_predict_train)

## get our estimates
data_predict_test_hat = forest.predict( data_summary_test )

### give us model output
print ( forest.oob_score_ )
print ( forest.feature_importances_ )

### error bars and plot 
plt.scatter(data_predict_test,data_predict_test_hat)
plt.plot([-4.5, -2], [-4.5, -2], 'k--')
plt.show()
plt.savefig('theta.png')

# Calculate the variance:
variance_unbiased = fci.random_forest_error( forest, data_predict_test, data_summary_test )

# Plot error bars for predicted MPG using unbiased variance
plt.errorbar( data_predict_test, data_predict_test_hat, yerr=np.sqrt(variance_unbiased), fmt='o')
plt.plot([-4.5, -2], [-4.5, -2], 'k--')
plt.show()
plt.savefig('theta.png')

### get our variance and estimators 
data_summary_test = np.concatenate( (data_summary_test, stats),axis=0 )
real_predict = forest.predict( data_summary_test )

variance_unbiased = fci.random_forest_error( forest, real_predict, data_summary_test ) 
print ( real_predict[-1] )
print ( variance_unbiased[-1] )
