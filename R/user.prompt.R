### Name: user.prompt
### Title: Scan user input
### Aliases: user.prompt
### Keywords: file

## The function is currently defined as
user.prompt <- function () {
 ANSWER <- readline("\nType 'y' to continue or 'n' to quit: ")
      
if (substr(ANSWER, 1, 1) == "n")
         {stop("Function terminated by user.")}
     
}



