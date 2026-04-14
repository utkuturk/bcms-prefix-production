PennController.ResetPrefix(null);

// Preload the ZIP files (will start with the top one, so give practice items first)
//Preload("https://raw.githubusercontent.com/ccuonzo/Italian-Latinate-Roots/refs/heads/main")
PreloadZip("https://raw.githubusercontent.com/ccuonzo/Serbian2/refs/heads/main/practice.zip", "https://raw.githubusercontent.com/ccuonzo/Serbian2/refs/heads/main/fillerNWs1.zip", "https://raw.githubusercontent.com/ccuonzo/Serbian2/refs/heads/main/fillerWs.zip", "https://raw.githubusercontent.com/ccuonzo/Serbian2/refs/heads/main/primes.zip", "https://raw.githubusercontent.com/ccuonzo/Serbian2/refs/heads/main/targets.zip");
//DebugOff();

// Runtime summary stats shown at the end
let expTargetCorrect = 0;
let expTargetTotal = 0;
let expStartTimestamp = null;

CheckPreloaded()
// Sequence of trials, based on their labels
PennController.Sequence( "consentform" ,
    "playfile" ,
    "instructions" ,
    randomize("practice") ,
    "start" ,
   "preload1" ,
  shuffle(randomize("1Y"), randomize("1N")) ,
    "break" ,
   "preload2",
shuffle(randomize("2Y"), randomize("2N")) ,
    "break" ,
   "preload3",
   shuffle(randomize("3Y"), randomize("3N")) ,
    "Questionaire",
    "send" ,
    "final"
    )
                          // Notice "send" before "final" --- refers to SendResults below
                          // Adjust amount of blocks + breaks in the sequence
    
PennController.DebugOff()    // Uncomment when ready to publish

// Start with the consent form
PennController( "consentform" ,
    newHtml("consentformserbian.html")
      .settings.checkboxWarning("Required")
      .settings.radioWarning("Required")
      .settings.inputWarning("Required")
      .print()
      .log()
    ,
    newButton("SLAŽEM SE")
        .center()
        .print()
        .wait(getHtml("consentformserbian.html").test.complete().failure(getHtml("consentformserbian.html").warn()))
)
.log( "ID" , PennController.GetURLParameter("id") )     // This will add the ID at the end of the results line
.setOption("hideProgressBar", true)                     // We don't show the progress bar for this trial

// Check that participants are able to play a "test" audio file    
PennController( "playfile" ,
    newHtml("testaudioserbian.html")
      .print()
    ,
    newButton("Kliknite ovde ako ste čuli zvučni test.")
        .center()
        .print()
        .wait()
)
.log( "ID" , PennController.GetURLParameter("id") )     
.setOption("hideProgressBar", true)  

// Give instructions
PennController( "instructions" ,
    newHtml("instructionsserbian.html")
      .print()
    ,
    newKey("press3","3").wait()
)
.log( "ID" , PennController.GetURLParameter("id") )     
.setOption("hideProgressBar", true)

// This creates a trial (labeled "preload"---see Sequence above) that only moves on when all the resources
// (ie all the audio files) used by the trials labeled ExpBlock1, etc. have been preloaded
PennController.CheckPreloaded("1Y" , "1N").label("preload1")
PennController.CheckPreloaded("2Y", "2N").label("preload2")
PennController.CheckPreloaded("3Y", "3N" ).label("preload3")
//PennController.CheckPreloaded().label("preloadall")
// --- PRACTICE ITEMS ---      
// This generates trials using the file Practiceitems.csv from chunk_includes
PennController.Template( "Practiceb2.csv" ,
    row => PennController( "practice" ,     // all these trials will be labeled 'practice' (see Sequence above)
        newText("F", "F: Prava reč")          // We create the Text elements in cache (not printed yet)
          .settings.bold()
        ,
        newText("J", "J: Besmislena reč")
          .settings.bold()
        ,
        newCanvas("text", 500, 100)    //
            .settings.center()
            .settings.add( 0 , 0 ,  getText("F") )
            .settings.add( 380 , 0 ,  getText("J") )    // aligned to the right edge
            .print()                                    // Prints the Canvas along with its Text elements
        ,
        newTimer("timer1",400-200*Math.random())    // This gives you a random ISI (here shorter for practice items since feedback stays for 900)
            .start()
            .wait()
        ,
        newAudio("prime", "https://raw.githubusercontent.com/ccuonzo/Serbian2/refs/heads/main/" + row.PrimeSoundfile)
            .settings.log("play","end")     // Logging when it starts and ends playing
            .play()
        ,
        newKey("answerPrime", "FJ")  // Respond by pressing F or J key
            .settings.log("all")
            .wait()
            .setVar("keyPressedprime")
            .settings.disable()
        ,
        newVar("keyPressedprime", "None")       
            .settings.global()
            .set( getKey("answerPrime") )  // This logs what the answer was
        ,
        getAudio("prime")                   // The key has been pressed (cf 'wait' on newKey above)
            .wait("first")                  // Now wait until audio has ended *if has not ended yet*
        ,
        getKey("answerPrime")
            .test.pressed(row.answerP)
            .success(newText("successP", "Tačno!").settings.color("green").settings.bold().settings.center().print())
            .failure(newText("failureP", "Netačno, molimo vas obratite više pažnje!").settings.color("red").settings.bold().settings.center().print())
        ,
        newTimer("timer2", 1500)
            .start()
            .wait()
        ,
        getText("successP").remove()
        ,
        getText("failureP").remove()
        ,
        newAudio("target", "https://raw.githubusercontent.com/ccuonzo/Serbian2/refs/heads/main/" + row.TargetSoundfile)
            .settings.log("play","end")
            .play()
        ,
        newKey("answerTarget", "FJ")  // Respond by pressing F or J key
            .settings.log("all")
            .wait()
            .setVar("keyPressedtarget")
            .settings.disable()
        ,
        newVar("keyPressedtarget", "None")       
             .settings.global()
             .set( getKey("answerTarget") )  // This logs what the answer was
        ,
        getAudio("target")
             .wait("first")
        ,
        getKey("answerTarget")
             .test.pressed(row.answerT)
             .success(newText("successT", "Tačno!").settings.color("green").settings.bold().settings.center().print())
             .failure(newText("failureT", "Netačno, molimo vas obratite više pažnje!").settings.color("red").settings.bold().settings.center().print())
        ,
        newTimer("timer3", 1500)
          .start()
          .wait()
    )
    .log("ID" , PennController.GetURLParameter("id") )
    .log("PrimeSoundfile"  , row.PrimeSoundfile  )
    .log("TargetSoundfile" , row.TargetSoundfile )
    .log("PCorrectAnswer", row.answerP)
    .log("TCorrectAnswer", row.answerT)
    .log("KeyPrime" , getVar("keyPressedprime") )   // key press prime  
    .log("KeyTarget" , getVar("keyPressedtarget") )  // key press target    
)

PennController( "start" ,
    newFunction("initSummaryStats", () => {
        expTargetCorrect = 0;
        expTargetTotal = 0;
        expStartTimestamp = Date.now();
    }).call()
    ,
    newHtml("startserbian.html")
      .print()
    ,
    newKey("press7","7").wait()
)
.setOption("hideProgressBar", true)

// --- EXPERIMENTAL ITEMS ---      
// This generates trials using the file Expitems.csv from chunk_includes
PennController.Template( "ExperimentalItemsb3.csv" ,
    row => PennController( row.Sequence ,           // These trials will be labeled from the BlockNum column; starting with the first block.
        newVar("TrialN", 0)
             .settings.global()
             .set(v => v+1 )
        ,
        newVar("ITI", 0)                // This will store the ITI
          .settings.global()            // Global so we can use it inside log below
          .set( v => Date.now() )       
        ,
        newText("F", "F: Prava reč")
          .settings.bold()
        ,
        newText("J", "J: Besmislena reč").settings.bold()
        ,
        newCanvas("text", 500, 100).settings.center()
            .settings.add( 0 , 0 ,  getText("F") )
            .settings.add( 380 , 0 ,  getText("J") )    // aligned to the right edge
            .print()                        // Prints the Canvas along with its Text elements
        ,
        newTimer("timer4", 600-200*Math.random())
            .start()
            .wait()
        ,
        getVar("ITI")                   // After newTimer
          .set( v => Date.now() - v )   // Set it to current timestamp - previous timestamp (v)
        ,
        newVar("primeRT", 0)                // This will store the RT for the prime
          .settings.global()                // Global so we can use it inside log below
          .set( v => Date.now() )           // Set it to the timestamp immediately before audio.play
        ,
        newAudio("prime", "https://raw.githubusercontent.com/ccuonzo/Serbian2/refs/heads/main/" + row.PrimeSoundfile)
            .settings.log("play","end")
            .play()
        ,     
        newKey("answerPrime", "FJ")  // Respond by pressing F or J key
            .settings.log("all")
            .wait()
            .setVar("keyPressedprime")
            .settings.disable()
        ,
        newVar("keyPressedprime", "None")       
            .settings.global()
            .set( getKey("answerPrime") )  // This logs what the answer given was
        ,
        getVar("primeRT")                   // The key has been pressed (cf wait on newKey above)
          .set( v => Date.now() - v )       // Set it to current timestamp - previous timestamp (v)
        ,
        getAudio("prime")
            .wait("first")
        ,
        newVar("ISI", 0)                // This will store the ISI = time between prime and target
            .settings.global()            // Global so we can use it inside log below
            .set( v => Date.now() )       
        ,    
        newTimer("timer5", 600-200*Math.random())
            .start()
            .wait()
        ,
        getVar("ISI")                   // After newTimer
          .set( v => Date.now() - v )   // Set it to current timestamp - previous timestamp (v)
        ,
        newVar("targetRT", 0)               // Same as above, but for target this time
          .settings.global()
          .set( v => Date.now() )
        ,
        newAudio("target", "https://raw.githubusercontent.com/ccuonzo/Serbian2/refs/heads/main/" + row.TargetSoundfile)
            .settings.log("play","end")
            .play()
        ,       
        newKey("answerTarget", "FJ")  // Respond by pressing F or J key
            .settings.log("all")
            .wait()
            .setVar("keyPressedtarget")
            .settings.disable()
        ,
        newVar("keyPressedtarget", "None")       
            .settings.global()
            .set( getKey("answerTarget") )  // This logs what the answer was
        ,
        newFunction("incrementTargetTotal", () => {
            expTargetTotal += 1;
        }).call()
        ,
        getKey("answerTarget")
            .test.pressed( String(row.tword).trim().toLowerCase() === "word" ? "F" : "J" )
            .success(
                newFunction("incrementTargetCorrect", () => {
                    expTargetCorrect += 1;
                }).call()
            )
        ,
        getVar("targetRT")
            .set( v => Date.now() - v )
        ,
        getAudio("target")
            .wait("first")
            )
    .log("ID" , PennController.GetURLParameter("id") )
    .log("PrimeSoundfile"  , row.PrimeSoundfile  )
    .log("TargetSoundfile" , row.TargetSoundfile )
    .log("primeRT" , getVar("primeRT") )           // Will append the values of primeRT and
    .log("targetRT" , getVar("targetRT") )     // targetRT to all of the trial's results lines
    .log("ITI" , getVar("ITI") )                   // ITI
    .log("ISI" , getVar("ISI") )
    .log("KeyPrime" , getVar("keyPressedprime") )   // gives the key press prime  
    .log("KeyTarget" , getVar("keyPressedtarget") )  // gives the key press target
    .log("Rownum"  , row.Rownum )
    .log("Sequence" , row.Sequence )
    .log("ItemIdent" , row.ItemIdent )
    .log("Condition" , row.Condition )
    .log("Group" , row.Group )
    .log("prime" , row.prime )
    .log("target" , row.target )
    .log("pword" , row.pword )
    .log("tword" , row.tword )
    .log("TargetExpectedKey" , String(row.tword).trim().toLowerCase() === "word" ? "F" : "J" )
    .log("ptype" , row.ptype )
    .log("TrialN", getVar("TrialN"))
)


// This creates a trial labeled "break"; we show it in between the blocks.
PennController( "break" ,
      newText("Sada možete napraviti kratku pauzu")
        .settings.bold()
        .print()
      ,
      newButton("Hajde da nastavimo")
        .center()
        .print()
        .wait()
)

// Questionaire at end
PennController( "Questionaire" ,
    newHtml("questionaireitHTML", "questionaireserbian.html")
      .print()
      .settings.radioWarning("Required")
      .settings.inputWarning("Required")
      .settings.log()// this logs the answers in the html file
    ,
    newVar("finalCode", "")
      .settings.global()
    ,
    newFunction("computeFinalCodeBeforeSend", () => {
        const alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
        const random2 = () =>
            alphabet.charAt(Math.floor(Math.random() * alphabet.length)) +
            alphabet.charAt(Math.floor(Math.random() * alphabet.length));

        const accuracyPct = expTargetTotal > 0 ? Math.round((expTargetCorrect / expTargetTotal) * 100) : 0;
        const accuracyStr = String(accuracyPct).padStart(2, "0");

        const elapsedMs = expStartTimestamp ? Date.now() - expStartTimestamp : 0;
        const totalMinutes = Math.floor(elapsedMs / 60000);
        const hh = String(Math.floor(totalMinutes / 60)).padStart(2, "0");
        const mm = String(totalMinutes % 60).padStart(2, "0");

        return random2() + "-" + accuracyStr + "-" + random2() + "-" + hh + "-" + mm;
    }).call().setVar("finalCode")
    , 
    newButton("Pošaljite odgovore i završite")
        .center()
        .print()
        .wait(getHtml("questionaireitHTML").test.complete().failure(getHtml("questionaireitHTML").warn()))
)
.log( "ID" , PennController.GetURLParameter("id") )     
.log( "FinalCode" , getVar("finalCode") )
.setOption("hideProgressBar", true)                     // We don't show the progress bar for this trial

    
// This is necessary to send the results *before* showing the final (everlasting) screen
PennController.SendResults("send")

// This creates the final screen
PennController( "final" ,
    newHtml("finalserbian.html")
      .print()
    ,
    newText("finalCodeDisplay", "")
        .settings.text( getVar("finalCode") )
        .settings.bold()
        .print()
    ,
    newTimer("finalTimer",1)
        .wait()                 // This will wait forever, because the Timer was never started
)
.setOption("countsForProgressBar", false)
