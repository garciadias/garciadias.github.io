const USER_PROFILE_URL = "https://api.github.com/users/";
var myUserName = "garciadias";

function fillTemplateWithUserData(data){
    var location = data.location || 'The Universe';
    $('#gitHubAvatar').attr("src", data.avatar_url);
    $('#gitHubUser').text(data.name + ' as @' + data.login);
    $('#gitHubLocation').text(location);
}

function fillTemplateWithDefaultData(XMLHttpRequest, textStatus, errorThrown){
    $('#gitHubAvatar').attr("src", 'assets/default-avatar.png');
    $('#gitHubUser').text('_Name_ _UserName_');
    $('#gitHubLocation').text('Frikilandia');
    console.log(textStatus);
    console.log(errorThrown);
    console.log(XMLHttpRequest.responseText);
}

$.ajax({
    type: 'GET',
    url: USER_PROFILE_URL + myUserName,
    timeout: 5000,
    success: fillTemplateWithUserData,
    error: fillTemplateWithDefaultData
});