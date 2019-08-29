# GenometricCorrelation project evaluating two genometric annotations 
#	correlation. Also provides some service IRanges-related procedures. 
# (c) 2010-2019 Alexander Favorov, Leslie Cope, Yulia Medvedeva, 
#              Loris Mularoni, Vsevolod Makeev, Sarah Wheelan.
#
# usage example



if(!require("GenometriCorr",quietly=TRUE)
	&& compareVersion( packageDescription("GenometriCorr")[["Version"]], "1.1.9" ) < 0 )
{
	install.packages('GenometriCorr',repos='http://genometricorr.sourceforge.net/R/',type='source')
	require("rtracklayer")
}

if(!require("tcltk",quietly=TRUE))
{
	stop("Your R instalaltion is without tcltk module;\nplease refer to the manual:\nhttp://cran.r-project.org/doc/manuals/R-admin.html#Tcl_002fTk")
}

#require(tcltk)      # Load the TclTk package

sink('genometricorrwizard.log')

window.to.raise='main'

my.temp.config.file=tempfile()

cat('Tmpfile:',my.temp.config.file,'\n')

#config<-new("GenometriCorrConfig",system.file("extdata","template-config.ini", package = "GenometriCorr"))
config<-new("GenometriCorrConfig")
config$data$query=""
config$data$reference=""

#config$data$query<-
#	system.file("extdata", "UCSCcpgis_hg19.bed", package = "GenometriCorr")

#config$data$reference<-
#	system.file("extdata", "UCSCrefseqgenes_hg19.bed", package = "GenometriCorr")

sink(my.temp.config.file)
print(config)
sink()

report.name='' #filename for report
plot.name=''
vizu.name=''
report.name.browsed='' #filename for report
plot.name.browsed=''
vizu.name.browsed=''

folder.set.by.config.operation=FALSE



main.dialog<-NULL
ofiles.dialog<-NULL
econfig.dialog<-NULL

go.button<-NULL #to be visible; it is OK button of main dialog
of.OK.button<-NULL #to be visible; it is OK button of output files dialog

#fileformats<-c('bed','gff','wig','bedGraph','bed15','gff1','gff2','gff3','bigwig')
fileformats<-c('bed','gff','wig','bedGraph','bed15','gff1','gff2','gff3')

action.tcl <- tclVar(0)   # tclVar() creates a Tcl variable

while(window.to.raise!='nothing,all done!')
{
	
	if(window.to.raise=='go')
	{
		window.to.raise='nothing,all done!'
		if(!is.null(main.dialog))tkdestroy(main.dialog)
		if(!is.null(ofiles.dialog))tkdestroy(ofiles.dialog)
		if(!is.null(econfig.dialog))tkdestroy(econfig.dialog)
		print(config)
		print(paste(c("Output:",report.name,plot.name,vizu.name),sep=' '))
		print('startted ....')
		while(sink.number()) sink() #unsink
		result<-GenometriCorr:::run.config(config)
		if (report.name!='')
		{
			sink(report.name)
			print(result)
			sink()
		}
		if(plot.name!='')
		{
			graphical.report(result,pdffile=plot.name)
		}
		if(vizu.name!='')
		{
			visualize(result,pdffile=plot.name)
		}
		next
	}

	if(window.to.raise=='exit' || window.to.raise=='destroy')
	{
		window.to.raise='nothing,all done!'
		if(!is.null(main.dialog))tkdestroy(main.dialog)
		if(!is.null(ofiles.dialog))tkdestroy(ofiles.dialog)
		if(!is.null(econfig.dialog))tkdestroy(econfig.dialog)
		unlink(my.temp.config.file)
		while(sink.number()) sink() #unsink
		next
	}

	if(window.to.raise=='ofiles.destroy' || window.to.raise=='ofiles.exit')
	{
		window.to.raise='main'
		tkdestroy(ofiles.dialog)
		ofiles.dialog<-NULL
		next
	}

	if(window.to.raise=='econfig.destroy' || window.to.raise=='econfig.exit')
	{
		print('econfig.destroy')
		window.to.raise='main'
		tkdestroy(econfig.dialog)
		econfig.dialog<-NULL
		next
	}

	if (window.to.raise=='main')
	{
		if(is.null(main.dialog))
		{
			#creating maindialog

			main.dialog <- tktoplevel(padx=3,pady=12)  # Create a new toplevel window
			#action.tcl <- tclVar(0)   # tclVar() creates a Tcl variable
			#main.dialog <- tktoplevel()  # Create a new toplevel window

			#buttonframe=tkframe(main.dialog, padding=message('"3 3 12 12"'))
			buttonframe=tkframe(main.dialog)

			load.button <- tkbutton(buttonframe, text = "Load configuration file",
					command = function() tclvalue(action.tcl) <- 'load.config')

			edit.button <- tkbutton(buttonframe, text = "Edit configuration file",
					command = function() tclvalue(action.tcl) <- 'edit.config')

			save.button <- tkbutton(buttonframe, text = "Save configuration file",
					command = function() tclvalue(action.tcl) <- 'save.config')

			cof.button <- tkbutton(buttonframe, text = "Choose output files",
					command = function() tclvalue(action.tcl) <- 'choose.output.files')

			go.button <- tkbutton(buttonframe, text = "Go!",
					command = function() tclvalue(action.tcl) <- 'go', state='disabled')

			#tkconfigure(go.button,state='normal')
			exit.button <- tkbutton(buttonframe, text = "Exit",
					command = function() tclvalue(action.tcl) <- 'exit')


			# Capture the event "Destroy" (e.g. Alt-F4 in Windows) and when this happens,
			# assign 'destroy' to action.tcl
			tkbind(main.dialog, "<Destroy>", function() tclvalue(action.tcl) <- 'destroy' )

			tktitle(main.dialog) <- "GenometriCorr Wizard: Configure the Correlation"  # Give the window a title

			tkgrid(buttonframe,row=0,column=0,sticky='nwes')

			tkgrid.columnconfigure(main.dialog,0,weight=1)
			tkgrid.rowconfigure(main.dialog,0,weight=1)

			tkgrid(load.button,column=0,row=0,sticky='we')
			tkgrid.configure(load.button,padx=5,pady=5)	

			tkgrid(edit.button,column=0,row=1,sticky='we')
			tkgrid.configure(edit.button,padx=5,pady=5)	

			tkgrid(save.button,column=0,row=2,sticky='we')
			tkgrid.configure(save.button,padx=5,pady=5)	

			tkgrid(cof.button,column=0,row=3,sticky='we')
			tkgrid.configure(cof.button,padx=5,pady=5)

			tkgrid(go.button,column=0,row=4,sticky='we')
			tkgrid.configure(go.button,padx=5,pady=5)	

			tkgrid(exit.button,column=0,row=5,sticky='we')
			tkgrid.configure(exit.button,padx=5,pady=5)	
		}
		if ( (report.name!='' && plot.name!='' && vizu.name!='') && 
				!is.null(config$data$query) && !is.null(config$data$reference) ) 
		{	
			print ('go.button+')
			tkconfigure(go.button,state='normal')
		}
		else
		{
			print ('go.button-')
			tkconfigure(go.button,state='disabled')
		}
		tkraise(main.dialog)
		tkfocus(main.dialog)         # Place the focus to our main tk window
		tkwait.variable(action.tcl)
		action.string <- as.character(tclvalue(action.tcl))   # Get and coerce content of a Tcl variable

		print(action.string)
		


		# Test the result
		if (action.string == 'load.config') 
		{
			f<-as.character(tclvalue(tkgetOpenFile()))
			if (f!='')
			{
				newconfig<-try(new('GenometriCorrConfig',f),silent=TRUE)
				if (inherits(newconfig,'GenometriCorrConfig'))
				{
					config<-newconfig
					sink(my.temp.config.file)
					print(config)
					sink()
					#we change path to the congig's path
					setwd(dirname(f))
					folder.set.by.config.operation=TRUE
				}
				else
				 tk_messageBox(type = "ok", message = 'Something went wrong.')
			}
			window.to.raise='main'
			next
		}
		if (action.string == 'edit.config') 
		{
			window.to.raise='edit.config'
			next
		}
		if (action.string == 'save.config')
		{
			f<-as.character(tclvalue(tkgetSaveFile()))
			print(c('Save:',f))
			if (f!='')
			{
			 	sink(f)
				print(config)
				sink()
				setwd(dirname(f))
				folder.set.by.config.operation=TRUE
			}
			window.to.raise='main'
			next
		}
		if (action.string == 'choose.output.files') 
		{
			window.to.raise='choose.output.files'
			#window.to.raise='main'
			next
		}
		if (action.string == 'go') 
		{
			window.to.raise='go'
			next
		}
		if (action.string == 'exit') window.to.raise='exit'
		if (action.string == 'destroy') window.to.raise='destroy'
		next
	}

	if (window.to.raise=='edit.config')
	{
		if (is.null(econfig.dialog))
		{
			econfig.dialog <- tktoplevel(padx=3,pady=12)  # Create a new toplevel window

			econfigframe=tkframe(econfig.dialog)


			#action.tcl <- tclVar(0)   # tclVar() creates a Tcl variable

			ec.OK.button <- tkbutton(econfigframe, text = "OK",
					command = function() tclvalue(action.tcl) <- 'econfig.ok')#,state='disabled')
			

			ec.Cancel.button <- tkbutton(econfigframe, text = "Cancel",
					command = function() tclvalue(action.tcl) <- 'econfig.cancel')


			browse.query.button <- tkbutton(econfigframe, text = "Browse query",
					command = function() tclvalue(action.tcl) <- 'econfig.browse.query')
			browse.reference.button <- tkbutton(econfigframe, text = "Browse reference",
					command = function() tclvalue(action.tcl) <- 'econfig.browse.reference')
			browse.mapping.button <- tkbutton(econfigframe, text = "Browse mapping",
					command = function() tclvalue(action.tcl) <- 'econfig.browse.mapping')

			#q.label <- tklabel(econfigframe, text="Query:")
			#r.label <- tklabel(econfigframe, text="Reference:")
			edit.label <- tklabel(econfigframe,text="Edit the configuration file:")
			files.label <- tklabel(econfigframe, text="Choose the input files by browsing:")
			#empty.label <- tklabel(econfigframe, text="   ")
			#config.edit.field<-tktext(econfigframe,bg="white",font="courier")
			config.scroll<-tkscrollbar(econfigframe, repeatinterval=5, command=function(...)tkyview(config.edit.field,...))
			config.edit.field<-tktext(econfigframe,bg="white",font="courier",yscrollcommand=function(...)tkset(config.scroll,...))

			sink(my.temp.config.file)
			print(config)
			sink()
			connection <- file(my.temp.config.file)
			Lines<-readLines(connection)
			close(connection)
			for (li in Lines)
				tkinsert(config.edit.field,'end',paste(li,"\n",sep=''))


			tktitle(econfig.dialog) <- "GenometriCorr Wizard: edit the configuration file"  # Give the window a title

			tkgrid(econfigframe,row=0,column=0,sticky='nwes')

			#tkgrid.columnconfigure(econfig.dialog,0,weight=1)
			#tkgrid.rowconfigure(econfig.dialog,0,weight=1)
			#kgrid.columnconfigure(econfig.dialog,2,pad=10)

			tkgrid(edit.label,row=0,column=0,sticky='we')
			tkgrid.configure(edit.label,padx=5,pady=5)

			tkgrid(files.label,row=0,column=2,columnspan=3,sticky='we')
			tkgrid.configure(files.label,padx=5,pady=5)
			
			tkgrid(config.edit.field,row=1,column=0,rowspan=5,sticky='we')
			tkgrid(config.scroll,row=1,column=1,rowspan=5,sticky='we')
			
			#tkgrid(empty.label,row=0,column=2,sticky='we')
			#tkgrid.configure(empty.label,padx=15,pady=5)

			#tkgrid(q.label,row=0,column=3,sticky='we')
			#tkgrid.configure(q.label,padx=5,pady=5)

			tkgrid(browse.query.button,row=1,column=3,columnspan=2,sticky='we')
			tkgrid.configure(browse.query.button,padx=5,pady=5,sticky='we')

			#tkgrid(r.label,row=1,column=3,sticky='we')
			#tkgrid.configure(r.label,padx=5,pady=5)

			tkgrid(browse.reference.button,row=2,column=3,columnspan=2,sticky='we')
			tkgrid.configure(browse.reference.button,padx=5,pady=5,sticky='we')

			#tkgrid(m.label,row=2,column=3,sticky='we')
			#tkgrid.configure(m.label,padx=5,pady=5)

			tkgrid(browse.mapping.button,row=3,column=3,columnspan=2,sticky='we')
			tkgrid.configure(browse.mapping.button,padx=5,pady=5,sticky='we')

			tkgrid(ec.Cancel.button, row=5, column=3,sticky='we')
			tkgrid.configure(ec.Cancel.button,padx=5,pady=5)

			tkgrid(ec.OK.button,row=5,column=4,sticky='we')
			tkgrid.configure(ec.OK.button,padx=5,pady=5)

			# Capture the event "Destroy" (e.g. Alt-F4 in Windows) and when this happens,
			# assign 2 to done.tcl
			tkbind(econfig.dialog, "<Destroy>", function() tclvalue(action.tcl) <- 'econfig.destroy' )
		}
		#if (
		#	as.character(tclvalue(report.tcl))!='' 
		#	|| 
		#	as.character(tclvalue(plot.tcl))!='' 
		#	|| 
		#	as.character(tclvalue(plot.tcl))!='') 
		#		tkconfigure(of.OK.button,state='active')
		econfig.signalled<-FALSE
		while(!econfig.signalled)
		{
			tkraise(econfig.dialog)
			tkfocus(econfig.dialog)         # Place the focus to our main tk window
			tkwait.variable(action.tcl)
			action.string <- as.character(tclvalue(action.tcl))   # Get and coerce content of a Tcl variable

			print(action.string)

			if (action.string == 'econfig.browse.query')
			{
				f<-as.character(tclvalue(tkgetOpenFile()))
				print(c('Browse q:',f))
				#save window content
				Lines<-strsplit(tclvalue(tkget(config.edit.field,'0.0','end')),'\n')
				sink(my.temp.config.file)
				for (li in Lines[[1]])
					cat(li,'\n',sep='')
				sink()
				newconfig<-try(new('GenometriCorrConfig',my.temp.config.file),silent=TRUE)
				if (inherits(newconfig,'GenometriCorrConfig'))
					config<-newconfig
				else
				 tk_messageBox(type = "ok", message = 'Something went wrong with the window content.')
				if (f!='')
				{
					#change query, replace file
					config$data$query=f
					sink(my.temp.config.file)
					print(config)
					sink()
					#clear window
					tkdelete(config.edit.field,'0.0','end')
					#re-read window
					connection <- file(my.temp.config.file)
					Lines<-readLines(connection)
					close(connection)
					for (li in Lines)
						tkinsert(config.edit.field,'end',paste(li,"\n",sep=''))
				}
				econfig.signalled=TRUE
				window.to.raise='edit.config'
				next
			}
			if (action.string == 'econfig.browse.reference')
			{
				f<-as.character(tclvalue(tkgetOpenFile()))
				print(c('Browse r:',f))
				#save window content
				Lines<-strsplit(tclvalue(tkget(config.edit.field,'0.0','end')),'\n')
				sink(my.temp.config.file)
				for (li in Lines[[1]])
					cat(li,'\n',sep='')
				sink()
				newconfig<-try(new('GenometriCorrConfig',my.temp.config.file),silent=TRUE)
				if (inherits(newconfig,'GenometriCorrConfig'))
					config<-newconfig
				else
				 tk_messageBox(type = "ok", message = 'Something went wrong with the window content.')
				if (f!='')
				{
					#change reference, replace file
					config$data$reference=f
					sink(my.temp.config.file)
					print(config)
					sink()
					#clear window
					tkdelete(config.edit.field,'0.0','end')
					#re-read window
					connection <- file(my.temp.config.file)
					Lines<-readLines(connection)
					close(connection)
					for (li in Lines)
						tkinsert(config.edit.field,'end',paste(li,"\n",sep=''))
				}
				econfig.signalled=TRUE
				window.to.raise='edit.config'
				next
			}
			if (action.string == 'econfig.browse.mapping')
			{
				f<-as.character(tclvalue(tkgetOpenFile()))
				print(c('Browse map:',f))
				#save window content
				Lines<-strsplit(tclvalue(tkget(config.edit.field,'0.0','end')),'\n')
				sink(my.temp.config.file)
				for (li in Lines[[1]])
					cat(li,'\n',sep='')
				sink()
				newconfig<-try(new('GenometriCorrConfig',my.temp.config.file),silent=TRUE)
				if (inherits(newconfig,'GenometriCorrConfig'))
					config<-newconfig
				else
				 tk_messageBox(type = "ok", message = 'Something went wrong with the window content.')
				if (f!='')
				{
					#change reference, replace file
					config$data$mapping=f
					sink(my.temp.config.file)
					print(config)
					sink()
					#clear window
					tkdelete(config.edit.field,'0.0','end')
					#re-read window
					connection <- file(my.temp.config.file)
					Lines<-readLines(connection)
					close(connection)
					for (li in Lines)
						tkinsert(config.edit.field,'end',paste(li,"\n",sep=''))
				}
				econfig.signalled=TRUE
				window.to.raise='edit.config'
				next
			}
			if (action.string == 'econfig.cancel') 
			{
				econfig.signalled=TRUE
				window.to.raise='econfig.exit'
				next
			}
			if (action.string == 'econfig.ok') 
			{
				Lines<-strsplit(tclvalue(tkget(config.edit.field,'0.0','end')),'\n')
				sink(my.temp.config.file)
				for (li in Lines[[1]])
					cat(li,'\n',sep='')
				sink()
				newconfig<-try(new('GenometriCorrConfig',my.temp.config.file),silent=TRUE)
				if (inherits(newconfig,'GenometriCorrConfig'))
					config<-newconfig
				else
				 tk_messageBox(type = "ok", message = 'Something went wrong with the window content.')
				econfig.signalled=TRUE
				window.to.raise='econfig.exit'
				next
			}
			if (action.string == 'econfig.destroy') 
			{
				econfig.signalled=TRUE
				window.to.raise='econfig.destroy'
				next
			}
			if (action.string == 'destroy') #global destroy 
			{
				econfig.signalled=TRUE
				window.to.raise='destroy'
				next
			}
			if (action.string == 'exit') #global destroy 
			{
				econfig.signalled=TRUE
				window.to.raise='exit'
				next
			}
		} #econfig.signalled
		next
	}

	if (window.to.raise=='choose.output.files')
	{
		if (is.null(ofiles.dialog))
		{
			ofiles.dialog <- tktoplevel(padx=3,pady=12)  # Create a new toplevel window

			ofilesframe=tkframe(ofiles.dialog)


			#action.tcl <- tclVar(0)   # tclVar() creates a Tcl variable
			report.tcl <- tclVar(report.name)
			plot.tcl <- tclVar(plot.name)
			vizu.tcl <- tclVar(vizu.name)

			of.OK.button <- tkbutton(ofilesframe, text = "OK",
					command = function() tclvalue(action.tcl) <- 'ofiles.ok')#,state='disabled')
			

			of.Cancel.button <- tkbutton(ofilesframe, text = "Cancel",
					command = function() tclvalue(action.tcl) <- 'ofiles.cancel')
			#set.report.button <- tkbutton(ofilesframe, text = "Set",
			#		command = function() tclvalue(action.tcl) <- 'ofiles.set.report')
			#set.plot.button <- tkbutton(ofilesframe, text = "Set",
			#		command = function() tclvalue(action.tcl) <- 'ofiles.set.plot')
			#set.vizu.button <- tkbutton(ofilesframe, text = "Set",
			#		command = function() tclvalue(action.tcl) <- 'ofiles.set.vizu')
			browse.report.button <- tkbutton(ofilesframe, text = "Browse",
					command = function() tclvalue(action.tcl) <- 'ofiles.browse.report')
			browse.plot.button <- tkbutton(ofilesframe, text = "Browse",
					command = function() tclvalue(action.tcl) <- 'ofiles.browse.plot')
			browse.vizu.button <- tkbutton(ofilesframe, text = "Browse",
					command = function() tclvalue(action.tcl) <- 'ofiles.browse.vizu')

			report.label <- tklabel(ofilesframe, text="Text report:")
			plot.label <- tklabel(ofilesframe, text="Plot pdf:")
			vizu.label <- tklabel(ofilesframe, text="Visualization pdf:")


			report.entry <- tkentry(ofilesframe, textvariable=report.tcl ,width=50,relief='sunken')
			plot.entry <- tkentry(ofilesframe, textvariable=plot.tcl ,width=50,relief='sunken')
			vizu.entry <- tkentry(ofilesframe, textvariable=vizu.tcl ,width=50,relief='sunken')

			tktitle(ofiles.dialog) <- "GenometriCorr Wizard: choose the output files"  # Give the window a title

			tkgrid(ofilesframe,row=0,column=0,sticky='nwes')

			tkgrid.columnconfigure(ofiles.dialog,0,weight=1)
			tkgrid.rowconfigure(ofiles.dialog,0,weight=1)

			tkgrid(report.label,row=0,column=0,sticky='we')
			tkgrid.configure(report.label,padx=5,pady=5)

			tkgrid(plot.label,row=1,column=0,sticky='we')
			tkgrid.configure(report.label,padx=5,pady=5)

			tkgrid(vizu.label,row=2,column=0,sticky='we')
			tkgrid.configure(report.label,padx=5,pady=5)

			tkgrid(report.entry,row=0,column=1,columnspan=6,sticky='we')
			tkgrid.configure(report.entry,padx=5,pady=5)

			tkgrid(plot.entry,row=1,column=1,columnspan=6,sticky='we')
			tkgrid.configure(report.entry,padx=5,pady=5)

			tkgrid(vizu.entry,row=2,column=1,columnspan=6,sticky='we')
			tkgrid.configure(report.label,padx=5,pady=5)

			#tkgrid(set.report.button,row=0,column=7,sticky='we')
			#tkgrid.configure(set.report.button,padx=0,pady=5,sticky='we')

			#tkgrid(set.plot.button,row=1,column=7,sticky='we')
			#tkgrid.configure(set.plot.button,padx=0,pady=5)

			#tkgrid(set.vizu.button,row=2,column=7,sticky='we')
			#tkgrid.configure(set.vizu.button,padx=0,pady=5)

			tkgrid(browse.report.button,row=0,column=8,sticky='we')
			tkgrid.configure(browse.report.button,padx=5,pady=5,sticky='we')

			tkgrid(browse.plot.button,row=1,column=8,sticky='we')
			tkgrid.configure(browse.plot.button,padx=5,pady=5)

			tkgrid(browse.vizu.button,row=2,column=8,sticky='we')
			tkgrid.configure(browse.vizu.button,padx=5,pady=5)

			tkgrid(of.Cancel.button, row=3, column=1,sticky='we')
			tkgrid.configure(of.Cancel.button,padx=5,pady=5)

			tkgrid(of.OK.button,row=3,column=5,sticky='we')
			tkgrid.configure(of.OK.button,padx=5,pady=5)

			# Capture the event "Destroy" (e.g. Alt-F4 in Windows) and when this happens,
			# assign 2 to done.tcl
			tkbind(ofiles.dialog, "<Destroy>", function() tclvalue(action.tcl) <- 'ofiles.destroy' )
		}
		#if (
		#	as.character(tclvalue(report.tcl))!='' 
		#	|| 
		#	as.character(tclvalue(plot.tcl))!='' 
		#	|| 
		#	as.character(tclvalue(plot.tcl))!='') 
		#		tkconfigure(of.OK.button,state='active')
		ofiles.signalled=FALSE
		while(!ofiles.signalled)
		{
			tkraise(ofiles.dialog)
			tkfocus(ofiles.dialog)         # Place the focus to our main tk window
			tkwait.variable(action.tcl)
			action.string <- as.character(tclvalue(action.tcl))   # Get and coerce content of a Tcl variable

			print(action.string)

			if (action.string == 'ofiles.browse.report')
			{
				f<-as.character(tclvalue(tkgetSaveFile()))
				print(c('Of browse report:',f))
				if (f!='')
				{
					tclvalue(report.tcl)<-f
					report.name.browsed<-f
				}

				#if(!folder.set.by.config.operation && f!="")
				#	setwd(dirname(f))
				
				ofiles.signalled=TRUE
				window.to.raise='choose.output.files'
				next
			}
			if (action.string == 'ofiles.browse.plot')
			{
				f<-as.character(tclvalue(tkgetSaveFile()))
				print(c('Of browse plot:',f))
				if (f!='')
				{
					tclvalue(plot.tcl)<-f
					plot.name.browsed<-f
				}

				#if(!folder.set.by.config.operation)
				#	setwd(dirname(f))

				ofiles.signalled=TRUE
				window.to.raise='choose.output.files'
				next
			}
			if (action.string == 'ofiles.browse.vizu')
			{
				f<-as.character(tclvalue(tkgetSaveFile()))
				print(c('Of browse vizu:',f))
				if (f!='')
				{
					tclvalue(vizu.tcl)<-f
					vizu.name.browsed<-f
				}

				#if(!folder.set.by.config.operation)
				#	setwd(dirname(f))

				ofiles.signalled=TRUE
				window.to.raise='choose.output.files'
				next
			}
			if (action.string == 'ofiles.cancel') 
			{
				ofiles.signalled=TRUE
				window.to.raise='ofiles.exit'
				next
			}
			if (action.string == 'ofiles.ok') 
			{
				report.name=as.character(tclvalue(report.tcl))
				if (report.name.browsed!=tclvalue(report.tcl)) 
					if (!(length(strsplit(basename(report.name),'\\.')[[1]]) > 1 && tail(strsplit(basename(report.name),'\\.')[[1]],1) == 'txt'))
						report.name=paste(report.name,'.txt',sep='')
				plot.name=as.character(tclvalue(plot.tcl))
				if (plot.name.browsed!=tclvalue(plot.tcl)) 
					if (!(length(strsplit(basename(plot.name),'\\.')[[1]]) > 1 && tail(strsplit(basename(plot.name),'\\.')[[1]],1) == 'pdf'))
						plot.name=paste(plot.name,'.pdf',sep='')
				vizu.name=as.character(tclvalue(vizu.tcl))
				if (vizu.name.browsed!=tclvalue(vizu.tcl)) 
					if (!(length(strsplit(basename(vizu.name),'\\.')[[1]]) > 1 && tail(strsplit(basename(vizu.name),'\\.')[[1]],1) == 'pdf'))
						vizu.name=paste(vizu.name,'.pdf',sep='')
				if (vizu.name==plot.name && vizu.name!="")
				{
					if (!(length(strsplit(basename(vizu.name),'\\.')[[1]]) > 1 && tail(strsplit(basename(vizu.name),'\\.')[[1]],1) == 'pdf'))
						vizu.name=paste(vizu.name,'.vizu',sep='')
					else #pdf on end
						vizu.name=paste(head(strsplit(basename(vizu.name),'\\.')[[1]],-1),'.vizu.pdf',sep='')
				}

				ofiles.signalled=TRUE
				window.to.raise='ofiles.exit'
				next
			}
			if (action.string == 'ofiles.destroy') 
			{
				ofiles.signalled=TRUE
				window.to.raise='ofiles.destroy'
				next
			}
			if (action.string == 'destroy') #global destroy 
			{
				ofiles.signalled=TRUE
				window.to.raise='destroy'
				next
			}
			if (action.string == 'exit') #global destroy 
			{
				ofiles.signalled=TRUE
				window.to.raise='exit'
				next
			}
		} #ofiles.signalled
		next
	}
	stop("Unknown action.\n")
	next
}
while(sink.number()) sink() #unsink again

