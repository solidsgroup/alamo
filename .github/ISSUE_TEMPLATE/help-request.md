---
name: Help Request
about: Please use this template to request help compiling or running Alamo.
title: 'Help requested: '
labels: help wanted
assignees: bsrunnels

---

Fill out this template if you are encountering issues with running Alamo. Please consult https://solidsgroup.github.io/alamo/docs/Questions.html for guidance on how to ask an effective question.

1. System configuration information
    - What is your OS (linux distro, mac os), machine type (virtual, laptop, desktop, hpc)
    - Were you able to clone successfuly?
    - Which alamo branch are you using?

2. Please attach your configure script. Start with the exact command you use to run your configure script:

    ```
     ./configure [args]
    ```

    and the entire output of the configure script.

    ```
    Please paste output here. Please ensure that this uses code block 
    formatting for readability (using triple back-ticks)
    ```

3. Please attach the output of your make command. (Preferably as a file since this can be quite large.)

4. Send run information.
    - If you are working on your own branch of Alamo, please ensure that the file is stored in that branch. Specify its location.
    - If you are using your own input file, please attach that input file and any dependencies here.
    - Send the exact command that you are using. For instance
        
        ``` mpirun -np 4 ./bin/alamo-2d-g++ myinputfile ```

5. Summarize expected results vs your current results. 

6. Summarize what you have tried.
