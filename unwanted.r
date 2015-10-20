cat("*** Now running regression models. Check", regressionsFilename, "for the output\n")
        
      ##regressionsCsvFilename= file.path(regressionsDirectory, paste("regressions", model, "csv", sep="."))
      ##stop()

      sink(regressionsFilename, append=TRUE)
      for ( j in 1:length(regressionVariables ) ) {
        ## clCounter=clusterCounter
        clCounter=1
        for ( level in levels(melted.roistats$cluster) ) {
          cat("####################################################################################################\n")
          cat("*** Seed: ", seed, "\n")
          cat("*** Running regression model for ", level, "on", regressionVariables[[j]]$name, "\n")
          
          ##csvLine=paste(level, group, phase, glt, regressionVariables[[i]]$name, sep=",")
          regressionFormula=as.formula(paste("value",  "~", regressionVariables[[j]]$variable, sep=" "))
          cat("*** Regression formula: ")
          print.formula(regressionFormula)
          mdl.df.all=   melted.roistats[melted.roistats$cluster %in% level, ]
          mdl.df.male=  melted.roistats[melted.roistats$Gender=="M"   & melted.roistats$cluster %in% level, ]
          mdl.df.female=melted.roistats[melted.roistats$Gender=="F" & melted.roistats$cluster %in% level, ]

          ct.all=cor.test(mdl.df.all[, regressionVariables[[j]]$variable], mdl.df.all$value, method="spearman")
          cat("### BOTH GENDERS COMBINED\n")
          print(ct.all)

          ct.male=cor.test(mdl.df.male[, regressionVariables[[j]]$variable], mdl.df.male$value, method="spearman")
          cat("### MALE ONLY\n")
          print(ct.male)

          ct.female=cor.test(mdl.df.female[, regressionVariables[[j]]$variable], mdl.df.female$value, method="spearman")
          cat("### FEMALE ONLY\n")
          print(ct.female)

          mdl.df.rlm=mdl.df.all[, c("value", regressionVariables[[j]]$variable, "Gender")]
          regression.formula=as.formula(paste("value", "~", regressionVariables[[j]]$variable, sep=" "))
          lm.model.nogender=lm(regression.formula, data=mdl.df.rlm)
          print(summary(lm.model.nogender))


          regression.formula=as.formula(paste("value", "~", regressionVariables[[j]]$variable, "+", "Gender", sep=" "))
          lm.model.withgender=lm(regression.formula, data=mdl.df.rlm)
          print(summary(lm.model.withgender))
          
          ##csvLine=sprintf("%s,%s,%s,%.2f,%.2f,%.5f,%.2f,%s",
          ##  seed, makeGraphTitle(level), regressionVariables[[j]]$name, ct$statistic, ifelse(is.null(ct$parameter), NA, ct$parameter),
          ##  ct$p.value, ct$estimate, make.significance.indications(ct$p.value))
          csvLine=sprintf("%s,%s,%s,%s,%s, %.2f,%.2f,%.5f,%.2f,%s",
            seed, "all", makeGraphTitle(level), paste(clusters[clCounter, c("CM RL", "CM AP", "CM IS")], collapse=","), regressionVariables[[j]]$name, ct.all$statistic,
            ifelse(is.null(ct.all$parameter), NA, ct.all$parameter),
            ct.all$p.value, ct.all$estimate, make.significance.indications(ct.all$p.value))
          cat(csvLine, "\n")
          push(mystack, csvLine)


          csvLine=sprintf("%s,%s,%s,%s,%s, %.2f,%.2f,%.5f,%.2f,%s",
            seed, "male", makeGraphTitle(level), paste(clusters[clCounter, c("CM RL", "CM AP", "CM IS")], collapse=","), regressionVariables[[j]]$name, ct.male$statistic,
            ifelse(is.null(ct.male$parameter), NA, ct.male$parameter),
            ct.male$p.value, ct.male$estimate, make.significance.indications(ct.male$p.value))
          cat(csvLine, "\n")
          push(mystack, csvLine)

          csvLine=sprintf("%s,%s,%s,%s,%s, %.2f,%.2f,%.5f,%.2f,%s",
            seed, "female", makeGraphTitle(level), paste(clusters[clCounter, c("CM RL", "CM AP", "CM IS")], collapse=","), regressionVariables[[j]]$name, ct.female$statistic,
            ifelse(is.null(ct.female$parameter), NA, ct.female$parameter),
            ct.female$p.value, ct.female$estimate, make.significance.indications(ct.female$p.value))
          cat(csvLine, "\n")
          push(mystack, csvLine)
          
          ##stop()
          
          if (ct.all$p.value < 0.05 ) {
            ## make scatter plots for the ROIs that are significant
            imageDirectory=file.path(group.results.dir, seed)
            if ( ! file.exists(imageDirectory) ) {
              dir.create(imageDirectory)
            }
            imageFilename=file.path(imageDirectory, sprintf("%s.fwhm%0.1f.%s.%s.and.%s.pdf", gsub(" +", ".", level),  usedFwhm, task, seed, regressionVariables[[j]]$variable))
            message(paste("*** Creating", imageFilename, "\n"))
            
            ylabel="RSFC"
            graphTitle=makeGraphTitle(level)
            ##sprintf("%s", toupper(seed))
            my.base.size=18

            graph=ggplot(mdl.df.all, aes(x=eval(parse(text=regressionVariables[[j]]$variable)), y=value)) +
              theme_bw(base_size =  my.base.size) +
                geom_point() +
                  geom_smooth(method="rlm", color="black") +
                  scale_fill_brewer(palette="Set1") +
                    labs(title = graphTitle, y=ylabel, x=regressionVariables[[j]]$name) +
                      theme(legend.position="none",
                            ##panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            axis.title.x = element_text(size=my.base.size, vjust=0),
                            axis.title.y = element_text(size=my.base.size, vjust=0.4, angle =  90),
                            plot.title=element_text(size=my.base.size*1.2, vjust=1))
            ggsave(imageFilename)
          } ## end of if (ct.all$p.value < 0.05 )
        
          clCounter=clCounter+1
          cat("clCounter is now", clCounter, "\n")
        } ## end of for ( level in levels(melted.roistats$cluster) )

        ##stop("Stopping")
      } ## end of for ( level in levels(roistats.summary$cluster) )
      
      sink()
      ##stop()
