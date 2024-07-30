for /R "C:\Program Files\R" %%I in (bin\Rscript.exe) do @if exist "%%I" "%%I" -e "shiny::runApp('app.R', launch.browser = TRUE)"
