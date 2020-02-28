/*jslint browser: true*/
/*global $*/
$(document).ready(function () {

    $('button').click(function () {
        var form_data = $('form').serializeArray()
        console.log(form_data)

        var group1 = form_data[0].value
        var group2 = form_data[1].value
        var email = form_data[3].value
        console.log(email)
        var genes = form_data[2].value
        console.log(genes)

        var json1 = JSON.stringify(group1);
        var json2 = JSON.stringify(group2);
        json_genes = JSON.stringify(genes);

        console.log(email)

        var emailPattern = /^[a-zA-Z0-9._-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,4}$/;
        var validemail = emailPattern.test(email)

        console.log(validemail)

        if (!validemail) {
            alert('Please enter a valid email address.');
            return false;
        }

        if (group1 == "") {
            alert(' You did not select any cells for group 1')
        }
        if (group2 == "") {
            alert(' You did not select any cells for group 2')
        }

        if (group1 != "" && group2 != "") {
            alert('GROUP 1: \n' + group1 + '\n GROUP 2: \n' + group2 + ' \n Genes: \n' + genes)

            $.post("/submit", {
                // "contentType": "application/json",
                'data1': json1,
                'data2': json2,
                'email': email,
                'genes': json_genes
            });
            location.reload();
        }
        ;
    });

});