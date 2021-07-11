Last Updated 07/06/2021

# Connecting to Minion with SSH

    To connect to the Minion, use the following command.

        ssh minit@IP.address

    The password is minit. The terminal should prompt you for the password.

    Notes:
        1.) You need to be connected to the same network as the Minion, which may
        require VPN.


# Copying Remote Files->Local Using SCP

    To copy files from the Minion to a local directory using scp, use the following 
    command. This is good for Unix-like operating systems. For Windows, you can use
    pscp.exe

        scp -r username@remoteHost:/remote/dir/file.txt /local/dir/

    To close the connection, use the following command.

        exit

    Notes:
        1.) The -r flag is recursive and will copy subfolders to the specified 
        directory as well.
        2.) You need to be connected to the same network as the Minion, which may
        require VPN.
        3.) If either directory has spaces in the path, enclose it in quotes.
        4.) To find the path of the data folder on the Minion, connect via SSH
        and go up directories until you find a subdirectory called data. The 
        data for a particular experiment will be in the subfolder
        /data/name_of_experiment/name_of_sample/

    Example:

        scp -r minit@IP.address:/remote/dir/file.txt "/my local/dir/"


# Updating the Minion Remotely
    
    If the Minion has issues updating locally, it can be updated remotely. To 
    remotely update the Minion, connect via SSH and run the following command.

        sudo apt update

    Reboot the Minion.

    Notes:
        1.) To connect via SSH, see "Connecting to Minion with SSH" above. 


# Counting Subsequence Matches