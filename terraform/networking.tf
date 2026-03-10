# Virtual Network
resource "azurerm_virtual_network" "vnet" {
  name                = "psuai894-vnet"
  address_space       = var.vnet_address_space
  location            = var.region
  resource_group_name = azurerm_resource_group.rg.name
  tags                = var.tags
}

# Subnet for Web App VNet Integration
resource "azurerm_subnet" "webapp_subnet" {
  name                 = "webapp-subnet"
  resource_group_name  = azurerm_resource_group.rg.name
  virtual_network_name = azurerm_virtual_network.vnet.name
  address_prefixes     = [var.webapp_subnet_prefix]

  delegation {
    name = "webapp-delegation"
    service_delegation {
      name    = "Microsoft.Web/serverFarms"
      actions = ["Microsoft.Network/virtualNetworks/subnets/action"]
    }
  }
}

# Subnet for Private Endpoint
resource "azurerm_subnet" "pe_subnet" {
  name                 = "pe-subnet"
  resource_group_name  = azurerm_resource_group.rg.name
  virtual_network_name = azurerm_virtual_network.vnet.name
  address_prefixes     = [var.pe_subnet_prefix]
}

# Private DNS Zone for Azure Files
resource "azurerm_private_dns_zone" "file_dns_zone" {
  name                = "privatelink.file.core.windows.net"
  resource_group_name = azurerm_resource_group.rg.name
  tags                = var.tags
}

# Link DNS Zone to VNet
resource "azurerm_private_dns_zone_virtual_network_link" "vnet_link" {
  name                  = "vnet-link"
  resource_group_name   = azurerm_resource_group.rg.name
  private_dns_zone_name = azurerm_private_dns_zone.file_dns_zone.name
  virtual_network_id    = azurerm_virtual_network.vnet.id
  tags                  = var.tags
}

# Private Endpoint for Storage Account
resource "azurerm_private_endpoint" "storage_pe" {
  name                = "storage-pe"
  location            = var.region
  resource_group_name = azurerm_resource_group.rg.name
  subnet_id           = azurerm_subnet.pe_subnet.id

  private_service_connection {
    name                           = "storage-connection"
    private_connection_resource_id = azurerm_storage_account.storage_account.id
    subresource_names              = ["file"]
    is_manual_connection           = false
  }

  private_dns_zone_group {
    name                 = "storage-dns-zone-group"
    private_dns_zone_ids = [azurerm_private_dns_zone.file_dns_zone.id]
  }

  tags = var.tags
}
