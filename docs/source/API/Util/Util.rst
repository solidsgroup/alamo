.. role:: cpp(code)
   :language: c++

Util
====

Basic Input/Output
------------------

- :code:`Util::Message`

    To generate a message:

    .. code-block:: cpp

        Util::Message(INFO,"This is a message. Here's a number: ", 8);

    This will generate the following output: note that the file name, line number (XXX)
    and function name are included in the message.

    .. code-block:: bash

        MESSAGE: src/path/to/file.cpp:XXX (function_name) This is a message. Here's a number: 8

    You can concatenate an arbitrary number of strings in one message, including :code:`Set::Vector`, etc.

    .. WARNING::

        You must include the :code:`INFO` argument as the first argument, otherwise it will
        not compile!

    You should always use :code:`Util::Message` instead of :code:`std::cout` so that Alamo
    can manage parallel IO.

- :code:`Util::Warning`

    Use this to warn the user any time you suspect dangerous behavior. 
    
    .. code-block:: cpp

        Util::Warning(INFO,"This is a message. I can print numbers like this ", 8);
    
    It is used the same way as :code:`Util::Message`.

- :code:`Util::Abort`

    Use this to safely and cleanly abort the simulation any time an error occurs.
    You can include an error message in the same way as with :code:`Util::Message`.

Documentation
-------------

.. doxygennamespace:: Util
   :project: alamo

