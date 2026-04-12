# Storage Account
resource "azurerm_storage_account" "storage_account" {
  name                          = var.storage_account_name
  resource_group_name           = azurerm_resource_group.rg.name
  location                      = var.region
  account_tier                  = "Standard"
  account_replication_type      = "LRS"
  public_network_access_enabled = false

  network_rules {
    default_action             = "Deny"
    virtual_network_subnet_ids = [azurerm_subnet.webapp_subnet.id]
    bypass                     = ["AzureServices"]
  }

  tags = var.tags
}

# Storage Share (Azure Files)
resource "azurerm_storage_share" "web_app_share" {
  name                 = var.app_service_name
  storage_account_name = azurerm_storage_account.storage_account.name
  quota                = 50
}
