<!-- @author Simo Endre (@simo_endre)
     Navier Stoke fluid simulation
-->

<!DOCTYPE HTML>
<head>

    <meta charset="utf-8">
    <meta http-equiv="content-type" content="text/html">
    <meta name="description" content="Navier Stoke fluid simulation">
    <meta name="author" content="Simo Endre">
    <meta name="keywords" content="fluid, generative, navier, stoke, pattern, noise, art, javascript, experiments"/>
    <meta name="viewport" content="width=device-width, initial-scale=1.0, maximum-scale=1.0"/>
    <link rel="shortcut icon" type="image/x-icon" href="assets/favicon.ico">

   <!-- <script type="text/javascript" src="http://platform.tumblr.com/v1/share.js"></script>-->

    <!--<link href='http://fonts.googleapis.com/css?family=Dosis' rel='stylesheet' type='text/css'>-->
    <link rel="stylesheet" href="css/fontface.css"></link>
    <link rel="stylesheet" href="css/style.css"></link>

    <title>HTML5 Fluid Simulation</title>

</head>
<body>

<header>
    <span id="wrap"></span>
    <span id="arrow"></span>
    <h1>Navier Stoke fluid simulation</h1>
    <span class='header-open'> click on <label id="label">+</label> for more info </span>

    <div class="vertical-line"></div>

</header> 
<div class="wrapper">   
    <div class="context"></div>	
</div>    
<div id="gui"></div>
<script>

    function getIEVersion() {
        var rv = -1; // Return value assumes failure.
        if (navigator.appName == 'Microsoft Internet Explorer') {
            var ua = navigator.userAgent;
            var re  = new RegExp("MSIE ([0-9]{1,}[\.0-9]{0,})");
            if (re.test(ua) != null)
                rv = parseFloat( RegExp.$1 );
        }
        return rv;
    }


    function checkVersion() {
        var ver = getIEVersion();

        if ( ver != -1 ) {
            if (ver <= 9.0) {
                document.body.innerHTML = "";
                document.body.innerHTML =
                        "<p class='no-canvas'>" +
                                "You need a <a href='http://www.google.com/chrome'>modern browser</a> to view this experiment." +
                        "</p>";
            }
        }
    }

    checkVersion();

</script>

<script src="script/header.js"></script>
<script src="libs/dat.gui.min.js"></script>
<script src="script/FS.js"></script>
<script src="script/FluidSolver.js"></script>
<script src="script/Point.js"></script>
<script src="script/Main.js"></script>

</body>
</html>