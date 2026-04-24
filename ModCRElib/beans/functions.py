"""
Shared utility functions used across ModCRE modules.

Includes lightweight file/FASTA parsers, queue and remote-execution helpers,
ortholog output parsing, and small workflow bookkeeping utilities.
"""

import os, sys, re
import gzip
import pwd
import subprocess
import time

#-------------#
# Parsers     #
#-------------#

def parse_file(file_name, gz=False):
    """
    Stream a text file line-by-line without trailing newlines.

    Args:
        file_name (str): Input path.
        gz (bool): If True, read as gzipped text.

    Yields:
        str: Next line with final ``\\n`` removed.

    Raises:
        ValueError: If ``file_name`` does not exist.
    """
    if os.path.exists(file_name):
        # Initialize #
        f = None
        # Open file handle #
        if gz: f = gzip.open(file_name, "rt")
        else: f = open(file_name, "rt")
        # For each line... #
        for line in f:
            yield line.strip("\n")
        f.close()
    else:
        raise ValueError("Could not open file %s" % file_name)

def parse_fasta_file(file_name, gz=False, clean=True):
    """
    Parse FASTA records as ``(identifier, sequence)`` tuples.

    Args:
        file_name (str): FASTA file path.
        gz (bool): If True, read as gzipped text.
        clean (bool): Replace non-word chars/digits by ``X`` in sequences.

    Yields:
        tuple[str, str]: Header text without ``>`` and the uppercase sequence.
    """
    # Initialize #
    identifier = ""
    sequence = ""
    # For each line... #
    for line in parse_file(file_name, gz):
        if line == "": continue
        if line[0] == ">":
            if sequence != "":
                if clean:
                    sequence = re.sub("\W|\d", "X", sequence)
                yield (identifier, sequence)
            m = re.search("^>(.+)", line)
            identifier = m.group(1)
            sequence = ""
        else:
            sequence += line.upper()
    if clean:
        sequence = re.sub("\W|\d", "X", sequence)

    yield (identifier, sequence)


def fileExist(file):
    """
    Check whether a path points to an existing regular file.

    Args:
        file (str | None): Candidate file path.

    Returns:
        bool: True if ``file`` is not ``None`` and exists as a regular file.
    """
    if file is not None:
        return os.path.exists(file) and os.path.isfile(file)
    else:
        return False


#-------------#
# Write       #
#-------------#

def write(file_name=None, content=None):
    """
    Append one line of text to file or stdout.

    Args:
        file_name (str | None): Output file path. If ``None``, prints to stdout.
        content (str): Text payload written as ``content + "\\n"``.

    Returns:
        None.

    Raises:
        ValueError: If file append fails.
    """
    if file_name is not None:
        try:
            f = open(file_name, "a")
            f.write("%s\n" % content)
        except:
            raise ValueError("Could not create file %s" % file_name)
    else:
        sys.stdout.write("%s\n" % content)

#-------------#
# Cluster     #
#-------------#

def submit_command_to_queue(command, queue=None, max_jobs_in_queue=None, queue_file=None, dummy_dir="/tmp", submit="qsub", qstat="qstat"):
    """
    Submit a shell command to a local cluster queue wrapper.

    Supports two modes: direct pipe to ``submit`` or generation of a temporary
    queue script from ``queue_file`` with the command appended.

    Args:
        command (str): Command line to execute in the queued job.
        queue (str | None): Queue/partition name; if ``None`` submit to default.
        max_jobs_in_queue (int | None): Optional throttle based on ``qstat``.
        queue_file (str | None): Optional template script prepended to submission.
        dummy_dir (str): Directory used to write temporary ``submit_*.sh`` files.
        submit (str): Submission executable (e.g. ``qsub``).
        qstat (str): Queue-status executable used for throttling.

    Returns:
        None.
    """
    import hashlib


    if max_jobs_in_queue is not None:
      try:
        while number_of_jobs_in_queue(qstat) >= max_jobs_in_queue: time.sleep(5)
      except Exception as e:
        print("Queue error to get the number of jobs (%s). Continue submitting jobs"%e)

    cwd = os.path.join(dummy_dir,"sh")
    if not os.path.exists(cwd): os.makedirs(cwd)
    script= os.path.join(cwd,"submit_"+hashlib.sha224(command.encode()).hexdigest()+".sh")
    if queue_file is not None:
      fd=open(script,"w")
      with open(queue_file,"r") as queue_standard:
        data=queue_standard.read()
        fd.write(data)
        fd.write("%s\n\n"%(command))
      fd.close()
      queue_standard.close()
      if queue is not None:
        try:
          process = subprocess.check_output([submit," -q ",queue,script])
        except Exception as e:
          print(("Failed execution %s"%e))
          try: 
            os.system("%s -q %s %s" % (submit, queue,script))
          except:
            print("Failed submission")
      else:
        try:
          process = subprocess.check_output([submit,script])
        except Exception as e:
          print(("Failed execution %s"%e))
          try: 
            os.system("%s %s"% (submit,script))
          except:
            print("Failed submission")

    else:
      if queue is not None:
        try:
          process = subprocess.check_output(["echo",command,"|",submit," -q ",queue])
        except Exception as e:
          print(("Failed execution %s"%e)) 
          try:
            os.system("echo \"%s\" | %s -q %s" % (command, submit, queue))
          except:
            print("Failed submission")
      else:
        try:
          process = subprocess.check_output(["echo",command,"|",submit])
        except Exception as e:
          print(("Failed execution %s"%e)) 
          try:
            os.system("echo \"%s\" | %s" % (submit,command))
          except:
            print("Failed submission")

def submit_command_to_server(command, config, dummy_dir="/tmp"):
    """
    Submit a command through SSH/SFTP to the configured remote cluster.

    Uses ``config`` sections ``Cluster`` and ``Paths`` to stage a generated
    script and execute the remote submit command.

    Args:
        command (str): Command to run remotely.
        config: ``configparser.ConfigParser`` with remote credentials/settings.
        dummy_dir (str): Local/remote scratch base used for submission scripts.

    Returns:
        str: Remote job id when parsing succeeds, otherwise an ``Error ...`` string.
    """
    import hashlib
    import paramiko

    #create client object
    ssh_client=paramiko.SSHClient()
    ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())

    #Collect data from configuration
    queue             = config.get("Cluster","cluster_queue")
    max_jobs_in_queue = config.get("Cluster","max_jobs_in_queue")
    submit            = config.get("Cluster","cluster_submit")
    qstat             = config.get("Cluster","cluster_qstat")
    host              = config.get("Cluster","server_host")
    user              = config.get("Cluster","server_user")
    passwd            = config.get("Cluster","server_passwd")
    remote_folder     = config.get("Cluster","server_directory")
    remote_outdir     = config.get("Cluster","server_output")
    scripts_path      = config.get("Paths","scripts_path")
    queue_file        = os.path.join(scripts_path,config.get("Cluster","command_queue"))

    #open connection 
    ssh_client.connect(hostname = host, username = user, password = passwd)
    
    #open a dummy folder in server
    remote_dummy_dir      = os.path.join(remote_outdir,dummy_dir)
    remote_cwd            = os.path.join(remote_dummy_dir,"sh")
    print(("REMOTE SUBMISSION: remote dummy directory ..."+remote_dummy_dir))
    print(("REMOTE SUBMISSION: remote command directory  ..."+remote_cwd))
    #because the address is the same
    os.system("mkdir %s"%remote_dummy_dir)
    os.system("chmod -R 777 %s"%remote_dummy_dir)
    os.system("mkdir %s"%remote_cwd)
    os.system("chmod -R 777 %s"%remote_cwd)
    #try also from the server
    stdin, stdout, stderr = ssh_client.exec_command("mkdir "+remote_dummy_dir)
    stdin, stdout, stderr = ssh_client.exec_command("chmod -R 777 "+remote_dummy_dir)
    stdin, stdout, stderr = ssh_client.exec_command("mkdir "+remote_cwd)
    stdin, stdout, stderr = ssh_client.exec_command("chmod -R 777 "+remote_cwd)
    stdin, stdout, stderr = ssh_client.exec_command("cd "+remote_cwd)

    if not os.path.exists(dummy_dir):
        os.makedirs(dummy_dir)
    os.system("chmod -R 777 "+dummy_dir)
    cwd = os.path.join(dummy_dir,"sh")
    if not os.path.exists(cwd): 
        os.makedirs(cwd)
    os.system("chmod -R 777 "+cwd)
    script= os.path.join(cwd,"submit_"+hashlib.sha224(command.encode()).hexdigest()+".sh")
    remote_script= os.path.join(remote_cwd,"submit_"+hashlib.sha224(command.encode()).hexdigest()+".sh")
    print(("SERVER QUEUE FILE: "+queue_file))
    print(("SERVER SUBMISSION: "+script))
    print(("REMOTE SUBMISSION: "+remote_script))
    try:
      if queue_file is not None:
       #Home files
       print(("REMOTE QUEUE: "+script))
       fd=open(script,"w")
       with open(queue_file,"r") as queue_standard:
        data=queue_standard.read()
        fd.write(data)
        fd.write("%s\n\n"%(command))
       fd.close()
       queue_standard.close()
       #Remote upload files
       print("REMOTE CONNECTING SSH CLIENT")
       ftp_client = ssh_client.open_sftp()
       print("REMOTE CONNECTING COPY QUEUE-SCRIPT ON CLIENT")
       ftp_client.put(script,remote_script)
       print("REMOTE CLOSE CONNECTION")
       ftp_client.close()
       #Create remote command
       if queue != "None":
          exec_command="%s -p %s %s" % (submit, queue,remote_script)
       else:
          exec_command="%s %s"% (submit,remote_script)
       print(("COMMAND = "+exec_command))
      else:
       if queue != "None":
          exec_command="echo \"%s\" | %s -p %s" % (command, submit, queue)
       else:
          exec_command="echo \"%s\" | %s" % (submit,command)
       print(("COMMAND = "+exec_command))

      #Run remote command
      print(("SEND JOB: %s "%exec_command))
      stdin, stdout, stderr = ssh_client.exec_command(exec_command)
      print("CHECKING ERRORS OF SERVER EXECUTION")
      for line in stderr:
         error=line
         print(line)
      print("CHECKING OUTPUT OF SERVER EXECUTION")
      for line in stdout:
        print(line)
        if line.startswith("Submitted"): job=line.split()[-1]

      try:
        return job
      except:
        raise error

    except Exception as e:

      print(("Error on server connection: "+str(e)))
      return "Error on server connection: "+str(e)


def execute_in_remote(command,parameters,config,output_dir,logfile,waiting=True):    
    """
    Create a remote wrapper script, submit it, and optionally wait for completion.

    Completion is detected by polling ``logfile`` for the command/job marker line.

    Args:
        command (str): Python script name (e.g. ``scanner.py``).
        parameters (str): CLI arguments passed to the script.
        config: ``ConfigParser`` with remote cluster settings.
        output_dir (str): Existing job directory, shared between local/remote FS.
        logfile (str): Log file path where DONE markers are appended.
        waiting (bool): If True, poll until the command marker appears.

    Returns:
        str: Job identifier or error message from remote submission helpers.
    """
        
    remote_scripts_path = os.path.join(config.get("Cluster","server_directory"),"scripts")
    remote_python       = os.path.join(config.get("Cluster","server_python"),"python")
    shell_file=os.path.join(output_dir,os.path.basename(output_dir)+"_"+command.rstrip(".py")+".sh")
    submit_file=open(os.path.join(output_dir,os.path.basename(output_dir)+"_"+command.rstrip(".py")+".sh"),"w")
    submit_file.write("echo '#START JOB %s  ' >> %s\n"%(os.path.basename(output_dir), logfile))
    submit_file.write("chmod 777 %s \n"%logfile)
    submit_file.write("%s %s %s \n"%(remote_python,os.path.join(remote_scripts_path,command),parameters))
    submit_file.write("echo '%s  %s   DONE' >> %s\n"%(command.rstrip(".py"), os.path.basename(output_dir), logfile))
    submit_file.write("echo '#FINISHED JOB  %s' >> %s\n"%( os.path.basename(output_dir), logfile))
    submit_file.write("chmod 777 %s \n"%logfile)
    submit_file.write("chmod -R 777 %s \n"%output_dir)
    submit_file.close()
    dummy_dir = "dummy_"+command.rstrip(".py")

    keep=[]
    if os.path.exists(logfile) and os.path.isfile(logfile):
        logread=open(logfile,"r")
        for line in logread:
            data=line.split()
            if command.rstrip(".py") == data[0]  and os.path.basename(output_dir) == data[1]: continue
            if command.rstrip(".py") in line and os.path.basename(output_dir) in line: continue
            keep.append(line.strip())
        logread.close()

    if len(keep)>0:
       log=open(logfile,"w")
       for line in keep:
           log.write("%s\n"%line)
       log.close()
       os.system("chmod 777 %s"%logfile)
    else:
       if os.path.exists(logfile) and os.path.isfile(logfile): os.remove(logfile)

    iterate = True
    #timing to write in output
    n_steps=0

    print(("\n\n... execute in remote .... (iterate %s , waiting %s , JOB_ID %s)"%(iterate,waiting,os.path.basename(output_dir))))
    
    #job = submit_command_to_server(" sh "+shell_file,config,dummy_dir)
    job = submit_command_to_server(" bash "+shell_file,config,dummy_dir)

    if str(job).startswith("Error"):
       return job

    print(("...starting job .... (iterate %s , waiting %s , JOB_ID %s , JOB_SUBMIT %s )"%(iterate,waiting,os.path.basename(output_dir),str(job))))

    if waiting:
     print(("...starting iteration.... (iterate %s , waiting %s , JOB_ID %s , JOB_SUBMIT %s )"%(iterate,waiting,os.path.basename(output_dir),str(job))))
     while(iterate):
      n_steps = n_steps + 1
      time.sleep(5)
      if n_steps>1.e+2:
        n_steps=0
        print(("...waiting to finish...(JOB_ID %s , LOG %s , JOB %s)"%(os.path.basename(output_dir),logfile,str(job))))
      if os.path.exists(logfile) and os.path.isfile(logfile):
          os.system("chmod 777 %s"%logfile)
          time.sleep(2)
          if n_steps==0: print(("...open log file...(JOB_ID %s , LOG %s )"%(os.path.basename(output_dir),logfile)))
          try:
            logread=open(logfile,"r")
            for line in logread:
                data=line.split()
                if n_steps==0: print(("...LOGFILE DATA: search command %s and JOB_ID %s"%(command.rstrip(".py"),os.path.basename(output_dir))))
                if n_steps==0: print(("...LOGFILE DATA: "+line))
                if command.rstrip(".py") == data[0]  and os.path.basename(output_dir) == data[1]: 
                    print(("...LOGFILE DATA FOUND: "+line))
                    iterate=False
                if command.rstrip(".py") in line and os.path.basename(output_dir) in line: 
                    print(("...LOGFILE DATA FOUND: "+line))
                    iterate=False
            logread.close()
          except Exception as e:
            print(("Error while reading %s: %s"%(logfile,e)))
            iterate=False
            continue
    return job

def submit_command_to_server3(command, config, dummy_dir="/tmp"):
    """
    Python3-oriented variant of ``submit_command_to_server``.

    It reads ``command_queue3`` from configuration and submits through the same
    SSH/SFTP flow as the Python2 variant.

    Args:
        command (str): Command to run remotely.
        config: ``ConfigParser`` with remote credentials/settings.
        dummy_dir (str): Local/remote scratch base used for submission scripts.

    Returns:
        str: Remote job id when parsing succeeds, otherwise an ``Error ...`` string.
    """
    import hashlib
    import paramiko

    #create client object
    ssh_client=paramiko.SSHClient()
    ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())

    #Collect data from configuration
    queue             = config.get("Cluster","cluster_queue")
    max_jobs_in_queue = config.get("Cluster","max_jobs_in_queue")
    submit            = config.get("Cluster","cluster_submit")
    qstat             = config.get("Cluster","cluster_qstat")
    host              = config.get("Cluster","server_host")
    user              = config.get("Cluster","server_user")
    passwd            = config.get("Cluster","server_passwd")
    remote_folder     = config.get("Cluster","server_directory")
    remote_outdir     = config.get("Cluster","server_output")
    scripts_path      = config.get("Paths","scripts_path")
    queue_file        = os.path.join(scripts_path,config.get("Cluster","command_queue3"))

    #open connection 
    ssh_client.connect(hostname = host, username = user, password = passwd)
    print("REMOTE SUBMISSION: USING PYTHON3 ") 
    #open a dummy folder in server
    remote_dummy_dir      = os.path.join(remote_outdir,dummy_dir)
    remote_cwd            = os.path.join(remote_dummy_dir,"sh")
    print(("REMOTE SUBMISSION: remote dummy directory ..."+remote_dummy_dir))
    print(("REMOTE SUBMISSION: remote command directory  ..."+remote_cwd))
    #because the address is the same
    os.system("mkdir %s"%remote_dummy_dir)
    os.system("chmod -R 777 %s"%remote_dummy_dir)
    os.system("mkdir %s"%remote_cwd)
    os.system("chmod -R 777 %s"%remote_cwd)
    #try also from the server
    stdin, stdout, stderr = ssh_client.exec_command("mkdir "+remote_dummy_dir)
    stdin, stdout, stderr = ssh_client.exec_command("chmod -R 777 "+remote_dummy_dir)
    stdin, stdout, stderr = ssh_client.exec_command("mkdir "+remote_cwd)
    stdin, stdout, stderr = ssh_client.exec_command("chmod -R 777 "+remote_cwd)
    stdin, stdout, stderr = ssh_client.exec_command("cd "+remote_cwd)

    if not os.path.exists(dummy_dir):
        os.makedirs(dummy_dir)
    os.system("chmod -R 777 "+dummy_dir)
    cwd = os.path.join(dummy_dir,"sh")
    if not os.path.exists(cwd): 
        os.makedirs(cwd)
    os.system("chmod -R 777 "+cwd)
    script= os.path.join(cwd,"submit_"+hashlib.sha224(command.encode()).hexdigest()+".sh")
    remote_script= os.path.join(remote_cwd,"submit_"+hashlib.sha224(command.encode()).hexdigest()+".sh")
    print(("SERVER SUBMISSION: "+script))
    print(("REMOTE SUBMISSION: "+remote_script))
    try:
      if queue_file is not None:
       #Home files
       print(("REMOTE QUEUE: "+script))
       fd=open(script,"w")
       with open(queue_file,"r") as queue_standard:
        data=queue_standard.read()
        fd.write(data)
        fd.write("%s\n\n"%(command))
       fd.close()
       queue_standard.close()
       #Remote upload files
       print("REMOTE CONNECTING SSH CLIENT")
       ftp_client = ssh_client.open_sftp()
       print("REMOTE CONNECTING COPY QUEUE-SCRIPT ON CLIENT")
       ftp_client.put(script,remote_script)
       print("REMOTE CLOSE CONNECTION")
       ftp_client.close()
       #Create remote command
       if queue != "None":
          exec_command="%s -p %s %s" % (submit, queue,remote_script)
       else:
          exec_command="%s %s"% (submit,remote_script)
       print(("COMMAND = "+exec_command))
      else:
       if queue != "None":
          exec_command="echo \"%s\" | %s -p %s" % (command, submit, queue)
       else:
          exec_command="echo \"%s\" | %s" % (submit,command)
       print(("COMMAND = "+exec_command))

      #Run remote command
      print(("SEND JOB: %s "%exec_command))
      stdin, stdout, stderr = ssh_client.exec_command(exec_command)
      print("CHECKING ERRORS OF SERVER EXECUTION")
      for line in stderr:
         error=line
         print(line)
      print("CHECKING OUTPUT OF SERVER EXECUTION")
      for line in stdout:
        print(line)
        if line.startswith("Submitted"): job=line.split()[-1]

      try:
         return job
      except:
         raise error

    except Exception as e:

      print(("Error on server connection: "+str(e)))

      return "Error on server connection: "+str(e)


def execute_in_remote3(command,parameters,config,output_dir,logfile,waiting=True):    
    """
    Python3-oriented variant of ``execute_in_remote``.

    Args:
        command (str): Python script name.
        parameters (str): CLI arguments passed to the script.
        config: ``ConfigParser`` with remote cluster settings.
        output_dir (str): Existing job directory, shared between local/remote FS.
        logfile (str): Log file path where DONE markers are appended.
        waiting (bool): If True, poll until the command marker appears.

    Returns:
        str: Job identifier or error message from remote submission helpers.
    """
        
    remote_scripts_path = os.path.join(config.get("Cluster","server_directory"),"scripts")
    remote_python       = os.path.join(config.get("Cluster","server_python3"),"python")
    shell_file=os.path.join(output_dir,os.path.basename(output_dir)+"_"+command.rstrip(".py")+".sh")
    submit_file=open(os.path.join(output_dir,os.path.basename(output_dir)+"_"+command.rstrip(".py")+".sh"),"w")
    submit_file.write("echo '#START JOB %s  ' >> %s\n"%(os.path.basename(output_dir), logfile))
    submit_file.write("chmod 777 %s \n"%logfile)
    submit_file.write("%s %s %s \n"%(remote_python,os.path.join(remote_scripts_path,command),parameters))
    submit_file.write("echo '%s  %s   DONE' >> %s\n"%(command.rstrip(".py"), os.path.basename(output_dir), logfile))
    submit_file.write("echo '#FINISHED JOB  %s' >> %s\n"%( os.path.basename(output_dir), logfile))
    submit_file.write("chmod 777 %s \n"%logfile)
    submit_file.write("chmod -R 777 %s \n"%output_dir)
    submit_file.close()
    dummy_dir = "dummy_"+command.rstrip(".py")

    keep=[]
    if os.path.exists(logfile) and os.path.isfile(logfile):
        logread=open(logfile,"r")
        for line in logread:
            data=line.split()
            if command.rstrip(".py") == data[0]  and os.path.basename(output_dir) == data[1]: continue
            if command.rstrip(".py") in line and os.path.basename(output_dir) in line: continue
            keep.append(line.strip())
        logread.close()

    if len(keep)>0:
       log=open(logfile,"w")
       for line in keep:
           log.write("%s\n"%line)
       log.close()
       os.system("chmod 777 %s"%logfile)
    else:
       if os.path.exists(logfile) and os.path.isfile(logfile): os.remove(logfile)

    iterate = True
    #timing to write in output
    n_steps=0

    print(("\n\n... execute in remote .... (iterate %s , waiting %s , JOB_ID %s)"%(iterate,waiting,os.path.basename(output_dir))))
    print(("\n\n... execute IN PYTHON3 .... (iterate %s , waiting %s , JOB_ID %s)"%(iterate,waiting,os.path.basename(output_dir)))) 
    #job = submit_command_to_server3(" sh "+shell_file,config,dummy_dir)
    job = submit_command_to_server3(" bash "+shell_file,config,dummy_dir)

    if str(job).startswith("Error"):
       return job

    print(("...starting job .... (iterate %s , waiting %s , JOB_ID %s , JOB_SUBMIT %s )"%(iterate,waiting,os.path.basename(output_dir),str(job))))

    if waiting:
     print(("...starting iteration.... (iterate %s , waiting %s , JOB_ID %s , JOB_SUBMIT %s )"%(iterate,waiting,os.path.basename(output_dir),str(job))))
     while(iterate):
      n_steps = n_steps + 1
      time.sleep(5)
      if n_steps>1.e+2:
        n_steps=0
        print(("...waiting to finish...(JOB_ID %s , LOG %s , JOB %s)"%(os.path.basename(output_dir),logfile,str(job))))
      if os.path.exists(logfile) and os.path.isfile(logfile):
          os.system("chmod 777 %s"%logfile)
          time.sleep(2)
          if n_steps==0: print(("...open log file...(JOB_ID %s , LOG %s )"%(os.path.basename(output_dir),logfile)))
          try:
            logread=open(logfile,"r")
            for line in logread:
                data=line.split()
                if n_steps==0: print(("...LOGFILE DATA: search command %s and JOB_ID %s"%(command.rstrip(".py"),os.path.basename(output_dir))))
                if n_steps==0: print(("...LOGFILE DATA: "+line))
                if command.rstrip(".py") == data[0]  and os.path.basename(output_dir) == data[1]: 
                    print(("...LOGFILE DATA FOUND: "+line))
                    iterate=False
                if command.rstrip(".py") in line and os.path.basename(output_dir) in line: 
                    print(("...LOGFILE DATA FOUND: "+line))
                    iterate=False
            logread.close()
          except Exception as e:
            print(("Error while reading %s: %s"%(logfile,e)))
            iterate=False
            continue
    return job



def number_of_jobs_in_queue(qstat="qstat"):
    """
    Count queued jobs for the current Unix user.

    Args:
        qstat (str): Queue-status executable.

    Returns:
        int: Number of lines in ``qstat -u <user>`` containing the username.
    """

    # Initialize #
    user_name = get_username()

    process = subprocess.check_output([qstat, "-u", user_name])

    return len([line for line in process.decode().split("\n") if user_name in line])

def get_username():
    """
    Return the current Unix account name.

    Returns:
        str: Username from ``pwd.getpwuid(os.getuid())``.
    """

    return pwd.getpwuid(os.getuid())[0]



#--------------------------------#
# Parse with configuration       #
#--------------------------------#


def parse_best_orthologs(input_file,pdb_dir,config,rank):
    """
    Parse scanner ortholog text output into ModCRE cluster dictionaries.

    Each ``>`` block describes one DNA interval/hit. Monomer and dimer entries
    are ranked by score and truncated according to ``max_orthologs`` unless
    ``rank`` is False.

    Args:
        input_file (str): Ortholog text file (``orthologs_with_best_templates.txt``).
        pdb_dir (str): PDB data directory containing ``families.txt``.
        config: ``ConfigParser`` with ``Parameters/max_orthologs``.
        rank (bool): If False, keep all ranked hits per interval.

    Returns:
        list[dict]: Interval dictionaries with proteins, families, monomer/dimer
        threading tuples, and metadata used by web JSON serializers.
    """

    orthologs_list = []

    # Get parameters from configuration
    max_orthologs = int(config.get("Parameters","max_orthologs"))
    #store families per pdb_chain
    families={}
    for line in parse_file(os.path.join(pdb_dir, "families.txt")):
            if line.startswith("#"): continue
            pdb_chain, family = line.split(";")
            families[pdb_chain] = family


    inp = open(input_file,"r")
    for line in inp:
        if line.startswith(">"): 
            #initialize
            ortholog_dict={}
            fragment_dict={}
            monomers     =[]
            dimers       =[]
            #store data
            data=line.strip().split("|")
            interval = data[1]
            start    = int(interval.split("-")[0])
            fragment_dict.setdefault("start",start)
            end      = int(interval.split("-")[1])
            fragment_dict.setdefault("end",end)
            pval     = float(data[2])
            fragment_dict.setdefault("pval",pval)
            hit_name = data[3]
            fragment_dict.setdefault("hit_name",hit_name)
            fragment_dict.setdefault("proteins",set())
            fragment_dict.setdefault("monomer",set())
            fragment_dict.setdefault("dimer",set())
            fragment_dict.setdefault("families",set())
            fragment_dict.setdefault("pdb_chains",set())
            continue
        if line.startswith("//"):
            #print "NEW ORTHOLOG LIST"
            if len(monomers)>0:
               #print monomers
               monomer_sorted = [monomer[1] for monomer in sorted(monomers,key=lambda x: x[0])]
               if not rank: max_orthologs=len(monomer_sorted)
               for thread in monomer_sorted[:min(max_orthologs,len(monomer_sorted))]:
                  uid,gene,thread_file,score,d_score = thread.split(";")
                  fragment_dict["proteins"].add((uid,gene))
                  fragment_dict["monomer"].add((thread_file,score,d_score))
                  pdb_Hchain = thread_file.split(".")[1]
                  pdb_chain  = pdb_Hchain[0:4] + "_" + pdb_Hchain[6:7]
                  fragment_dict["pdb_chains"].add(pdb_chain)
                  if pdb_chain in families: fragment_dict["families"].add(families[pdb_chain])
            if len(dimers)>0:
               #print dimers
               dimer_sorted   = [dimer[1] for dimer in sorted(dimers,key=lambda x: x[0])]
               if not rank: max_orthologs=len(dimer_sorted)
               for thread in dimer_sorted[:min(max_orthologs,len(dimer_sorted))]:
                  #read thread files 
                  monomer_A,monomer_B = thread
                  uid_A,gene_A,thread_A,score_A,d_score_A = monomer_A.split(";")
                  uid_B,gene_B,thread_B,score_B,d_score_B = monomer_B.split(";")
                  #define fragment_dict for uidA
                  fragment_dict["proteins"].add((uid_A,gene_A))
                  fragment_dict["monomer"].add((thread_A,score_A,d_score_A))
                  pdb_Hchain = thread_A.split(".")[1]
                  pdb_chain  = pdb_Hchain[0:4] + "_" + pdb_Hchain[6:7]
                  fragment_dict["pdb_chains"].add(pdb_chain)
                  if pdb_chain in families: fragment_dict["families"].add(families[pdb_chain])
                  #define fragment_dict for uidB
                  fragment_dict["proteins"].add((uid_B,gene_B))
                  fragment_dict["monomer"].add((thread_B,score_B,d_score_B))
                  pdb_Hchain = thread_B.split(".")[1]
                  pdb_chain  = pdb_Hchain[0:4] + "_" + pdb_Hchain[6:7]
                  fragment_dict["pdb_chains"].add(pdb_chain)
                  if pdb_chain in families: fragment_dict["families"].add(families[pdb_chain])
                  #define fragment_dict for dimer
                  dimer=((thread_A,score_A,d_score_A),(thread_B,score_B,d_score_B))
                  idimer=((thread_B,score_B,d_score_B),(thread_A,score_A,d_score_A))
                  if dimer in fragment_dict["dimer"]: continue
                  if idimer in fragment_dict["dimer"]: continue
                  fragment_dict["dimer"].add(dimer)
            if len(dimers)>0 or len(monomers)>0:
                  fragment_dict["dimer"]      = [x for x in fragment_dict["dimer"] ]
                  fragment_dict["monomer"]    = [x for x in fragment_dict["monomer"] ]  
                  fragment_dict["proteins"]   = [x for x in fragment_dict["proteins"] ]
                  fragment_dict["pdb_chains"] = [x for x in fragment_dict["pdb_chains"] ]
                  fragment_dict["families"]   = [x for x in fragment_dict["families"] ]
                  ortholog_dict = fragment_dict
            if ortholog_dict != {}: orthologs_list.append(ortholog_dict)
            continue
        if not line.startswith("//") and not line.startswith(">"):
            monomer_A,monomer_B = line.strip().split()
            if monomer_A == "None" or len(monomer_A.split(";"))<5: monomer_A=None
            if monomer_B == "None" or len(monomer_B.split(";"))<5: monomer_B=None
            if monomer_A is not None:
                uid_A,gene_A,thread_A,score_A,d_score_A = monomer_A.split(";")
                if score_A != 0 and d_score_A != 0 : monomers.append((float(score_A),monomer_A))
            if monomer_B is not None:
                uid_B,gene_B,thread_B,score_B,d_score_B = monomer_B.split(";")
                if score_B != 0 and d_score_B != 0 : monomers.append((float(score_B),monomer_B))
            if monomer_B is not None and monomer_A is not None:
                if score_B != 0 and d_score_B != 0 and score_A != 0 and d_score_A != 0 :
                   dimers.append((min(float(score_A),float(score_B)),(monomer_A,monomer_B)))

    inp.close()

    return orthologs_list





#-------------#
# Iteration   #
#-------------#


def done_jobs(info_file):
    """
    Read finished job identifiers from an info/log file.

    Args:
        info_file (str): Path where first token per non-comment line is a job id.

    Returns:
        set[str]: Completed job ids.
    """
    done=set()
    if not fileExist(info_file): return done
    info=open(info_file,"r")
    for line in info:
        if line.startswith("#"):continue
        done.add(line.strip().split()[0])
    info.close()
    return done

def check_done(done,set_of_jobs):
    """
    Check whether there are pending jobs not present in the done set.

    Args:
        done (iterable[str]): Completed job paths or names.
        set_of_jobs (iterable[str]): Target job paths or names.

    Returns:
        bool: True if at least one target job is not done.
    """
    if len(set_of_jobs)<=0:return False
    done_set=set()
    for data in done:
        done_set.add(os.path.basename(data))
    testing_set=set()
    for data in set_of_jobs:
        testing_set.add(os.path.basename(data))
    if len(done)>0: return ( testing_set!= testing_set.intersection(done_set))
    else: return True


def check_submitted(submitted,pdb_files):
    """
    Check whether there are PDB jobs still not present in submitted records.

    Args:
        submitted (iterable[str]): Submitted job paths or names.
        pdb_files (iterable[str]): Expected PDB job paths or names.

    Returns:
        bool: True if at least one PDB file has not been submitted.
    """
    submitted_set=set()
    for data in submitted:
        submitted_set.add(os.path.basename(data))
    pdb_set=set()
    for data in pdb_files:
        pdb_set.add(os.path.basename(data))
    if len(submitted)>0: return (pdb_set!= pdb_set.intersection(submitted_set))
    else: return True


#-------------#
# Conversion  #
#-------------#



def reverse_dna(seq):
    """
    Compute reverse complement for canonical DNA alphabet ``A,C,G,T``.

    Args:
        seq (str): Input DNA sequence.

    Returns:
        str: Reverse-complement sequence.
    """
    reverse={"A":"T","T":"A","C":"G","G":"C"}
    rev=""
    for n in range(len(seq)-1,-1,-1):
        rev=rev+reverse[seq[n]]
    return rev


