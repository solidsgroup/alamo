.. role:: cpp(code)
   :language: c++

.. _simba:

======================================
:fas:`question-circle;fa-fw` Questions
======================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:



Alamo is a research code that is maintained by a small band of developers consisting of graduate students, a PI, and our `collaborators <https://github.com/solidsgroup/alamo/graphs/contributors>`_.
While we are happy to provide help and support current (and future) collaborators, we have limited resources and cannot provide comprehensive support.
Therefore, we ask that all contributors stick to these guidelines when requesting help.

Where to ask for help
=====================

#. Submit an `issue on GitHub <https://github.com/solidsgroup/alamo/issues>`_.
   This is the preferred method as it allows us to track issues in a systematic way.

#. Send a message on Slack. If you are a developer and would like access to the :code:`alamodevelopers` please let us know by email.

#. Send an email to `brunnels@iastate.edu <mailto:brunnels@iastate.edu>`_.
   


Submitting a help request
=========================

If you are running into issues with the code, please use the following guidelines when submitting your help request.
Remember to send all the information needed in order to exactly reproduce the error that you are getting.

:bdg-success-line:`Every help request should come with complete, comprehensive, step-by-step instructions on how to exactly reproduce your issue.`

1.  Send all the information that is unique to your system.
    This includes the following:

   * :bdg-primary:`System configuration`
     What version of linux are you using, is this a laptop, virtual machine, etc.
     Have you installed all the necessary :ref:`dependencies<Installing dependencies>`?

   * :bdg-primary:`Alamo clone information`
     Were you able to clone Alamo successfully?
     If not: please send the *exact command* that you used to try to clone Alamo, along with the *complete, exact output* that you received.

   * :bdg-primary:`Alamo branch`
     Make sure that your local copy of Alamo is consistent with a branch on GitHub.

2. Send the complete **configure script output**: this includes

   * :bdg-primary:`./configure command` - this is the exact command that you used to run the configure script, e.g.

     .. code::

        ./configure --dim=2 --debug

   * :bdg-primary:`./configure output` - this is the **entire output** of your configure command, e.g.

     .. code::

        Python              3.8
        Current Branch      development
        Current Hash        b09f26ee5661fcd47d0f76e6ca00f4b5c1e276c4
        Dimension           3
        Debug Mode          False
        OpenMP              False
        ...
        ...

     Terminal screenshots are OK but searchable text is preferred.
     The configure script is designed to include a lot of diagnostic information.

     :bdg-danger-line:`ALWAYS` submit your configure script output!
     
3. Send your :bdg-primary:`make command` and the exact :bdg-primary:`output of your make command`.

   Note that you can pipe the output of your make command to a file (myfile.txt) by

   .. code::

      make > myfile.txt

   or

   .. code::

      make | tee myfile.txt
                      
   if you want to have the output go both to your terminal and to the file.
   (This works for any terminal output.)

4. Send run information.
   This includes:

   * The exact :bdg-primary:`input file` that you are using.
     If the input file is stored in the repository (this is preferred), then just make sure that the repository is up to 
     Otherwise, you should send your exact input file along with the help request.

   * The exact :bdg-primary:`command` that you are using.
     For instance,

     .. code::

        mpirun -np 2 ./bin/alamo-2d-g++ myinputfile

5. Summarize :bdg-primary:`expected results` vs :bdg-primary:`your current results`.
   Are you just trying to get alamo to compile?
   Are you trying to develop a new capability?
   Are you trying to create a new regression test?
   Are you trying to reproduce a result from a paper?
   The more context you can give - even if you are a current member - will expedite the request for help.

6. Summarize :bdg-primary:`what you have tried`.
   Do not expect anyone to help you if you have not made an honest effort to solve the problem yourself.
   You should have a complete summary of everything you have tried before requesting personalized help.
   

