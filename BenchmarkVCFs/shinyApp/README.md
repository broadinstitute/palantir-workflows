
Shiny App for viewing BenchmarkVcf resulting csv files
====


This app can be used to view the data contained in the csv file that is the result of running the BenchmarkVcf script. It can be used to two modes:

1. local: 
In R studio, open the `app.R` file and click on "Run App" button in the top right of the editor window (if you haven't installed the shiny package you may have to install it using `install.packages("shiny")`)
2. server:
You may create a user in shiny.io and publish this app there by clicking on the publish app button (right next to the Run App button) and following the instructions. (The publishing step will complain that After the app is published, anyone will be able to use it and view their results using this resource (access may be limited depending on your subscription). A version of this app is currently published [here](https://yfarjoun.shinyapps.io/shiny_vis_app/). 

In either case, the first thing you'll need to do is to upload the csv file: Click "Browse" and select the csv that you downloaded from the BenchmarkVcf run. 
Now you can look around. if you click a point in a plot or a row in a table your browser should open a gcloud console page pointing at a xml file. this is an IGV session that contains the vcfs that were used to arrive at the datum you clicked. 
