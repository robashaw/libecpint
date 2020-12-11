
.. _program_listing_file__Users_robertshaw_devfiles_libecpint_src_makelist.py:

Program Listing for File makelist.py
====================================

|exhale_lsh| :ref:`Return to documentation for file <file__Users_robertshaw_devfiles_libecpint_src_makelist.py>` (``/Users/robertshaw/devfiles/libecpint/src/makelist.py``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: py

   import sys
   max_am = int(sys.argv[1])
   prefix = str(sys.argv[2])
   
   file = open(prefix + "/qlist.txt", "w")
   file.write(prefix + "/generated/ecpint_gen.cpp\n")
   for j in range(max_am+1):
       for i in range(j+1):
           for k in range(max_am+1):
               if j == i == k == max_am:
                   file.write(prefix + "/generated/Q" + str(i) + str(j) + str(k) + ".cpp") 
               else:
                   file.write(prefix + "/generated/Q" + str(i) + str(j) + str(k) + ".cpp\n")
   file.close()
