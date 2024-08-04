
# Connecting to Computerome Server

I used the Computerome server to take advantage of its powerful computing resources. It provides more processing power and storage than my local system, which is essential for handling large datasets and running complex analyses efficiently.
here are 3 ways to login to Computerome's HPC environment: 
1. SSH via terminal
1. SSH via client
1. Virtual desktop
   
The following explanation is based on SSH via Terminal as a login service. 

1. On a Mac or Linux machine open a new Terminal and execute the following command (replace `username` with the one you received from Computerome). 
   ```
   $ ssh username@ssh.computerome.dk
   ```
1. Next, enter **the temporary password** you received via SMS. 
1. Now you will receive an Entrust Identity mobile notification to confirm. 
1. After you have confirmed the notification you will be logged into the HPC environment.

The disk quota of my personal home directory (`/home/people/user`) was 10 GB. But luckily Computerome provided me with another directory with much larger disk quota. Therefore, I changed the current working directory using this command:
```
$cd /home/projects/cu_10160/people/elehos
