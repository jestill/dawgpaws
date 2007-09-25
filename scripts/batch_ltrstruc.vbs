''----------------------------------------------------+
'' batch_ltrstruc.vbs                                 |
''----------------------------------------------------+
''  AUTHOR: James C. Estill                           |
'' CONTACT: JamesEstill_at_gmail.com                  |
'' STARTED: 12/13/2006                                |
'' DESCRIPTION:                                       |
''  I think that this will automate the running       |
''  of the LTR_Struc program. This prevents the       |
''  user from needing to answer the questions         |
''  and hit the enter key for each BAC.               |
''                                                    |
'' From the command line this can be launched         |
'' with the Cscript program                           |
'' Cscript.exe //NoLogo batch_ltrstruc.vbs | cmd.exe  |
''----------------------------------------------------+
Wscript.Stdout.WriteLine "LTR_STRUC_1_1"
Wscript.Sleep 1000
Wscript.Stdout.WriteLine "Y"
Wscript.Sleep 1000
'' This select N for using the standard settings
Wscript.Stdout.WriteLine "N"
'' This selects 1 for the most stringent
Wscript.Stdout.WriteLine "1"