**Monitoring script**

The Monitoring script is a script that collects memory, cpu and disk information about a task while it's running.
If you control your cromwell configuration, you can simply add it to your cromwell config as documented [here](https://cromwell.readthedocs.io/en/stable/wf_options/Google/).
If you cannot change the cromwell config (for example, if you are running on Terra), you can add monitoring to your task by adding 

    File monitoring_script
    
to the `inputs` of your task, and 
    
    bash ~{monitoring_script} > monitoring.log &
    
to the top of the `command` section of your task, and 

    File monitoring_log = "monitoring.log"
    
to the `output` section of your task. 

You can either hard code the monitoring script location, or provide it as input to your task.

A working monitoring script is provided here as `cromwell_monitoring_script.sh` it is designed to work 
on almost all linux docker images.

When your task is done you will have access to the monitoring.log file. 
You can plot the evolution of the various metrics with the supplied R script `plot_monitoring_data.R` with the command

    ./plot_monitoring_data monitoring.log monitoring.pdf
    
Which will take monitoring.log as input and generate the monitoring.pdf as output.

An example can be seen in the examples directory.
 