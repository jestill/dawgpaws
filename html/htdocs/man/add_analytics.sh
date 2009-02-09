# Add the google analytics script to html docs
# usage cat my_file.html | add_analytics.sh > my_file.html
#sed 's/<\/body>/soul/'
sed -i 's/<\/body>/<script src="http:\/\/www.google-analytics.com\/urchin.js"\n type="text\/javascript">" \n<\/script>\n<script type="text\/javascript">\n_uacct = "UA-2207628-5";\nurchinTracker();\n<\/script>\n<\/body>/' $1
#
# NOTES:
# Need to append the following
#
#<script src="http://www.google-analytics.com/urchin.js"
# type="text/javascript">
#</script>
#type="text/javascript">"
#</script>
#<script type="text/javascript">
#_uacct = "UA-2207628-5";
#urchinTracker();
#</script>

